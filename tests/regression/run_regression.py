#!/usr/bin/env python3
"""Run the ILLUMINA kernel on a synthetic case and compare the results.

The kernel reads ``illumina.in`` from the current directory when it gets
no command line arguments, so the binary is executed with no arguments
and the case directory as working directory (the plain, single-pointing
code path). The six summary numbers are parsed from the
``<basenm>.out`` file (after the last ``====`` line) and the contribution
map ``<basenm>_pcl.bin`` is parsed at full float32 precision. Both are
compared against ``reference.json`` in the case directory.

Exit status is 0 when every value is within tolerance, 1 otherwise.
Use ``--update`` to (re)write the reference from the current binary.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

LABELS = [
    "irdirect",   # Direct irradiance from sources (W/m**2/nm)
    "irrdirect",  # Direct irradiance from reflexion (W/m**2/nm)
    "direct",     # Direct radiance from sources (W/str/m**2/nm)
    "rdirect",    # Direct radiance from reflexion (W/str/m**2/nm)
    "cloud",      # Cloud radiance (W/str/m**2/nm)
    "diffuse",    # Diffuse radiance (W/str/m**2/nm)
]


def read_basenm(case):
    with open(os.path.join(case, "illumina.in")) as f:
        f.readline()
        return f.readline().split()[0]


def parse_out(path):
    """Return the six numbers that follow the last '====' line."""
    with open(path) as f:
        lines = [ln.strip() for ln in f]
    idx = max(i for i, ln in enumerate(lines) if "====" in ln)
    values = []
    for ln in lines[idx + 1:]:
        if not ln:
            continue
        try:
            values.append(float(ln))
        except ValueError:
            pass  # a label line
    if len(values) != len(LABELS):
        raise RuntimeError(
            "expected %d numbers after the last ==== line in %s, found %d"
            % (len(LABELS), path, len(values))
        )
    return dict(zip(LABELS, values))


def parse_los(path):
    """Return the line-of-sight step count and the furthest horizontal distance.

    A DUG comparison found the step count a cleaner observable of terrain
    blocking than any radiance: where the observer horizon was mishandled,
    the sight line ran out to the sources instead of stopping at the ridge,
    and the step count showed it plainly while the radiance only showed a
    ratio. Returns (None, None) when the log has no progression lines, which
    is the case for a quiet run.
    """
    steps, furthest = 0, None
    with open(path, errors="replace") as f:
        for line in f:
            if "Progression along the line of sight" in line:
                steps += 1
            elif "Horizontal dist. line of sight" in line:
                # "  Horizontal dist. line of sight =   29.5993938      m"
                try:
                    furthest = float(line.rsplit("=", 1)[1].replace("m", "").strip())
                except (IndexError, ValueError):
                    pass
    if steps == 0:
        return None, None
    return steps, furthest


def read_bin(path):
    """Read an ILLUMINA 2-D binary; return (nbx, nby, dict '(i,j)'->value).

    Keys use 1-based Fortran indices (i along x, j along y). Only non-zero
    values are returned; the map is sparse (one value per lamp pixel).
    """
    raw = np.fromfile(path, dtype=np.uint32)
    nbx, nby = int(raw[1]), int(raw[2])
    vals = np.frombuffer(raw[4:].tobytes(), dtype=np.float32).reshape(-1, 3)[:, 1]
    if vals.size != nbx * nby:
        raise RuntimeError("bad record count in %s" % path)
    out = {}
    for k in np.nonzero(vals)[0]:
        j = nby - int(k) // nbx
        i = int(k) % nbx + 1
        out["%d,%d" % (i, j)] = float(vals[k])
    return nbx, nby, out


def run_case(binary, case, timeout):
    binary = os.path.abspath(binary)
    basenm = read_basenm(case)
    for ext in (".out", "_pcl.bin"):
        try:
            os.remove(os.path.join(case, basenm + ext))
        except FileNotFoundError:
            pass
    t0 = time.time()
    with open(os.path.join(case, "run.log"), "w") as log:
        proc = subprocess.run([binary], cwd=case, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
    wall = time.time() - t0
    if proc.returncode != 0:
        sys.exit("kernel exited with status %d (see %s/run.log)" % (proc.returncode, case))
    result = {
        "summary": parse_out(os.path.join(case, basenm + ".out")),
        "wall_time_s": round(wall, 3),
    }
    los_steps, los_dist = parse_los(os.path.join(case, "run.log"))
    if los_steps is not None:
        result["los_steps"] = los_steps
        if los_dist is not None:
            result["los_max_dist_m"] = round(los_dist, 2)
    pcl = os.path.join(case, basenm + "_pcl.bin")
    if os.path.exists(pcl):
        nbx, nby, nz = read_bin(pcl)
        result["pcl_shape"] = [nbx, nby]
        result["pcl_nonzero"] = nz
    return result


def close(a, b, rtol, atol):
    return abs(a - b) <= atol + rtol * abs(b)


def compare(result, ref, rtol, atol):
    rows = []
    ok = True
    for k in LABELS:
        r, n = ref["summary"][k], result["summary"][k]
        good = close(n, r, rtol, atol)
        ok &= good
        rows.append((k, r, n, good))
    # Only when this run produced a log. run_angles.py drives the kernel
    # itself and compares the result through here, so it has no log of its
    # own; the step count is checked on the runs this script makes.
    if "los_steps" in ref and result.get("los_steps") is not None:
        r, n = ref["los_steps"], result["los_steps"]
        good = (n == r)
        ok &= good
        rows.append(("los_steps", r, n, good))
    ref_pcl = ref.get("pcl_nonzero")
    if ref_pcl is not None:
        new_pcl = result.get("pcl_nonzero", {})
        for k in sorted(set(ref_pcl) | set(new_pcl), key=lambda s: tuple(map(int, s.split(",")))):
            r, n = ref_pcl.get(k, 0.0), new_pcl.get(k, 0.0)
            good = close(n, r, rtol, atol)
            ok &= good
            rows.append(("pcl(%s)" % k, r, n, good))
        if ref.get("pcl_shape") != result.get("pcl_shape"):
            ok = False
            rows.append(("pcl_shape", str(ref.get("pcl_shape")), str(result.get("pcl_shape")), False))
    return ok, rows


def fmt(v):
    return "%14.7e" % v if isinstance(v, float) else "%14s" % v


def print_table(rows, rtol):
    print("%-16s %14s %14s   %s" % ("quantity", "reference", "current", "status"))
    for name, r, n, good in rows:
        rel = ""
        if isinstance(r, float) and isinstance(n, float) and r != 0:
            rel = "  rel=%.1e" % (abs(n - r) / abs(r))
        print("%-16s %s %s   %s%s" % (name, fmt(r), fmt(n), "ok" if good else "MISMATCH", rel))
    print("tolerance: rtol=%g" % rtol)


def ensure_case(case):
    if os.path.exists(os.path.join(case, "illumina.in")):
        return
    args = []
    ref_path = os.path.join(case, "reference.json")
    if os.path.exists(ref_path):
        with open(ref_path) as f:
            args = json.load(f).get("make_case_args", [])
    print("Case inputs missing; generating with make_case.py", *args)
    subprocess.check_call([sys.executable, os.path.join(HERE, "make_case.py"), "--out", case] + args)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", default=os.path.join(REPO, "bin", "illumina"))
    ap.add_argument("--case", default=os.path.join(HERE, "case_small"))
    ap.add_argument("--update", action="store_true", help="rewrite reference.json from this run")
    ap.add_argument("--rtol", type=float, default=1e-5)
    ap.add_argument("--atol", type=float, default=0.0)
    ap.add_argument("--timeout", type=float, default=1800.0)
    ap.add_argument(
        "--make-case-args",
        default=None,
        help="arguments recorded in reference.json for regenerating the case (with --update)",
    )
    args = ap.parse_args(argv)

    case = os.path.abspath(args.case)
    ensure_case(case)
    result = run_case(args.binary, case, args.timeout)
    ref_path = os.path.join(case, "reference.json")
    print("binary: %s   wall time: %.2f s" % (args.binary, result["wall_time_s"]))

    if args.update:
        ref = dict(result)
        ref["binary"] = os.path.abspath(args.binary)
        ref["created"] = time.strftime("%Y-%m-%d %H:%M:%S")
        try:
            ref["git_commit"] = subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"], cwd=REPO, text=True
            ).strip()
        except Exception:
            pass
        old_args = []
        if os.path.exists(ref_path):
            with open(ref_path) as f:
                old_args = json.load(f).get("make_case_args", [])
        ref["make_case_args"] = args.make_case_args.split() if args.make_case_args else old_args
        with open(ref_path, "w") as f:
            json.dump(ref, f, indent=2, sort_keys=True)
            f.write("\n")
        print("Reference written to", ref_path)
        for k in LABELS:
            print("  %-10s %.3e" % (k, result["summary"][k]))
        return 0

    if not os.path.exists(ref_path):
        sys.exit("no reference.json in %s; run with --update first" % case)
    with open(ref_path) as f:
        ref = json.load(f)
    ok, rows = compare(result, ref, args.rtol, args.atol)
    print_table(rows, args.rtol)
    print("reference wall time: %.2f s (%s)" % (ref.get("wall_time_s", float("nan")), ref.get("created", "?")))
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
