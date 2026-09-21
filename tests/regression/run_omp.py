#!/usr/bin/env python3
"""Check the OpenMP build of the ILLUMINA kernel against the serial build.

The script copies the inputs of a case into temporary work directories
and runs:

1. the serial binary (``bin/illumina``),
2. the OpenMP binary (``bin/illumina_omp``) with ``OMP_NUM_THREADS=1``:
   ``<basenm>_result.txt`` and ``<basenm>_pcl.bin`` must be byte-identical
   to the serial run,
3. the OpenMP binary with every other thread count of ``--threads``
   (default 4): the six results of ``<basenm>_result.txt`` and the
   non-zero pixels of ``<basenm>_pcl.bin`` must agree with the serial run
   within ``--rtol`` (default 1e-4). The summation order of the OpenMP
   reduction differs, so bit-identity is not expected here.

The wall time of every run and the largest relative difference of every
threaded run are printed. Exit status is 0 on PASS and 1 on FAIL.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)

import run_regression as rr  # noqa: E402
from run_angles import copy_inputs, read_bytes  # noqa: E402

RESULT_KEYS = [
    "direct_irradiance_sources",
    "direct_irradiance_reflection",
    "direct_radiance_sources",
    "direct_radiance_reflection",
    "cloud_radiance",
    "diffuse_radiance",
]


def parse_result(path):
    """Return the key=value lines of a ``_result.txt`` file as floats."""
    out = {}
    with open(path) as f:
        for ln in f:
            if "=" in ln:
                k, v = ln.split("=", 1)
                out[k.strip()] = float(v)
    missing = [k for k in RESULT_KEYS if k not in out]
    if missing:
        raise RuntimeError("%s lacks %s" % (path, ", ".join(missing)))
    return out


def run(binary, cwd, threads, log, timeout):
    env = dict(os.environ)
    if threads is not None:
        env["OMP_NUM_THREADS"] = str(threads)
    t0 = time.time()
    with open(log, "w") as fh:
        proc = subprocess.run(
            [os.path.abspath(binary)], cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env, timeout=timeout
        )
    wall = time.time() - t0
    if proc.returncode != 0:
        sys.exit("kernel exited with status %d (see %s)" % (proc.returncode, log))
    return wall


def rel_diff(a, b):
    if b == 0.0:
        return 0.0 if a == 0.0 else float("inf")
    return abs(a - b) / abs(b)


def compare_tol(ref_dir, new_dir, basenm, rtol, atol):
    """Compare the six results and the pcl pixels; return (ok, maxrel, rows)."""
    ok = True
    maxrel = 0.0
    rows = []
    ref = parse_result(os.path.join(ref_dir, basenm + "_result.txt"))
    new = parse_result(os.path.join(new_dir, basenm + "_result.txt"))
    for k in RESULT_KEYS:
        good = rr.close(new[k], ref[k], rtol, atol)
        ok &= good
        maxrel = max(maxrel, rel_diff(new[k], ref[k]))
        rows.append((k, ref[k], new[k], good))
    ref_pcl = rr.read_bin(os.path.join(ref_dir, basenm + "_pcl.bin"))[2]
    new_pcl = rr.read_bin(os.path.join(new_dir, basenm + "_pcl.bin"))[2]
    for k in sorted(set(ref_pcl) | set(new_pcl), key=lambda s: tuple(map(int, s.split(",")))):
        r, n = ref_pcl.get(k, 0.0), new_pcl.get(k, 0.0)
        good = rr.close(n, r, rtol, atol)
        ok &= good
        maxrel = max(maxrel, rel_diff(n, r))
        rows.append(("pcl(%s)" % k, r, n, good))
    return ok, maxrel, rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--serial", default=os.path.join(REPO, "bin", "illumina"), help="serial binary")
    ap.add_argument("--omp", default=os.path.join(REPO, "bin", "illumina_omp"), help="OpenMP binary")
    ap.add_argument("--case", default=os.path.join(HERE, "case_small"))
    ap.add_argument(
        "--threads",
        default="1,4",
        help="comma separated thread counts (default 1,4); 1 is checked for byte identity, the others within --rtol",
    )
    ap.add_argument("--rtol", type=float, default=1e-4)
    ap.add_argument("--atol", type=float, default=0.0)
    ap.add_argument("--timeout", type=float, default=7200.0)
    ap.add_argument("--keep", action="store_true", help="keep the temporary work directories")
    ap.add_argument("--verbose", action="store_true", help="print every compared value")
    args = ap.parse_args(argv)

    case = os.path.abspath(args.case)
    rr.ensure_case(case)
    basenm = rr.read_basenm(case)
    threads = [int(t) for t in args.threads.split(",") if t.strip()]
    if 1 not in threads:
        threads.insert(0, 1)

    work = tempfile.mkdtemp(prefix="illum_omp_")
    ok = True
    try:
        ser_dir = os.path.join(work, "serial")
        copy_inputs(case, ser_dir)
        t_ser = run(args.serial, ser_dir, None, os.path.join(work, "serial.log"), args.timeout)
        print("serial            %-40s wall %8.2f s" % (args.serial, t_ser))

        for nt in threads:
            omp_dir = os.path.join(work, "omp%d" % nt)
            copy_inputs(case, omp_dir)
            t_omp = run(args.omp, omp_dir, nt, os.path.join(work, "omp%d.log" % nt), args.timeout)
            speedup = t_ser / t_omp if t_omp > 0 else float("inf")
            line = "OMP_NUM_THREADS=%-2d %-40s wall %8.2f s  speedup %5.2f" % (nt, args.omp, t_omp, speedup)
            if nt == 1:
                same = True
                for suffix in ("_result.txt", "_pcl.bin"):
                    a = read_bytes(os.path.join(ser_dir, basenm + suffix))
                    b = read_bytes(os.path.join(omp_dir, basenm + suffix))
                    if a != b:
                        same = False
                        print("  MISMATCH: %s differs between serial and 1 thread" % (basenm + suffix))
                ok &= same
                print(line + ("  byte-identical" if same else "  NOT IDENTICAL"))
            else:
                good, maxrel, rows = compare_tol(ser_dir, omp_dir, basenm, args.rtol, args.atol)
                ok &= good
                print(line + "  max rel diff %.2e  %s" % (maxrel, "ok" if good else "MISMATCH"))
                if args.verbose or not good:
                    rr.print_table(rows, args.rtol)
    finally:
        if args.keep:
            print("work directories kept in", work)
        else:
            shutil.rmtree(work, ignore_errors=True)

    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
