#!/usr/bin/env python3
"""Check the multi-pointing loop of the ILLUMINA kernel.

The kernel accepts an optional third argument, an angles list file with
one ``elevation_deg azimuth_deg`` pointing per line. It loads the domain
once and loops over the pointings. This script checks that the loop
gives the same results as separate single-pointing processes:

1. three single runs, each in a temporary copy of the case with line 17
   of ``illumina.in`` rewritten (``30 45`` = the reference pointing,
   ``10 120``, ``60 300``), binary called with no arguments;
2. one looped run with an angles file that lists the three pointings
   (plus comment and blank lines), binary called as
   ``illumina illumina.in <basenm>.out angles.txt``;
3. every looped ``_result.txt`` and ``_pcl.bin`` must be byte-identical
   to the matching single run;
4. the first looped pointing must match ``reference.json``;
5. the combined ``<basenm>_results.txt`` must equal the three single
   ``_result.txt`` files joined with one blank line;
6. the wall times of the three single runs and of the looped run are
   printed;
7. one more looped run with the fourth argument ``nomaps`` must give
   byte-identical ``_result.txt`` files and no ``_pcl.bin`` file.

Exit status is 0 on PASS and 1 on FAIL.
"""
import argparse
import json
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

POINTINGS = [(30.0, 45.0), (10.0, 120.0), (60.0, 300.0)]
ANGLE_LINE = 17  # 1-based line of 'Elevation viewing angle ; Azimuthal viewing angle'


def angle_tag(value):
    """File name tag of an angle: one decimal, '.'->'p', leading '-'->'m'."""
    s = "%.1f" % value
    if s.startswith("-"):
        s = "m" + s[1:]
    return s.replace(".", "p")


def pointing_root(basenm, elev, azim):
    return "%s_e%s_a%s" % (basenm, angle_tag(elev), angle_tag(azim))


def is_input(name):
    if name in ("reference.json",) or name.endswith((".out", ".log")):
        return False
    if name.endswith("_pcl.bin") or "_result" in name:
        return False
    return True


def copy_inputs(case, dest):
    os.makedirs(dest)
    for name in os.listdir(case):
        src = os.path.join(case, name)
        if os.path.isfile(src) and is_input(name):
            shutil.copy(src, dest)


def rewrite_angle_line(parfile, elev, azim):
    with open(parfile) as f:
        lines = f.readlines()
    lines[ANGLE_LINE - 1] = (
        "%-30s ! Elevation viewing angle ; Azimuthal viewing angle\n" % ("%g %g" % (elev, azim))
    )
    with open(parfile, "w") as f:
        f.writelines(lines)


def run(binary, cwd, args, log, timeout):
    t0 = time.time()
    with open(log, "w") as fh:
        proc = subprocess.run([binary] + args, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, timeout=timeout)
    wall = time.time() - t0
    if proc.returncode != 0:
        sys.exit("kernel exited with status %d (see %s)" % (proc.returncode, log))
    return wall


def read_bytes(path):
    with open(path, "rb") as f:
        return f.read()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", default=os.path.join(REPO, "bin", "illumina"))
    ap.add_argument("--case", default=os.path.join(HERE, "case_small"))
    ap.add_argument("--rtol", type=float, default=1e-5)
    ap.add_argument("--atol", type=float, default=0.0)
    ap.add_argument("--timeout", type=float, default=1800.0)
    ap.add_argument("--keep", action="store_true", help="keep the temporary work directories")
    args = ap.parse_args(argv)

    binary = os.path.abspath(args.binary)
    case = os.path.abspath(args.case)
    rr.ensure_case(case)
    basenm = rr.read_basenm(case)
    with open(os.path.join(case, "reference.json")) as f:
        ref = json.load(f)

    work = tempfile.mkdtemp(prefix="illumina_angles_")
    ok = True
    try:
        # 1. three single-pointing runs
        single = []
        t_single = 0.0
        for k, (elev, azim) in enumerate(POINTINGS):
            d = os.path.join(work, "single_%d" % k)
            copy_inputs(case, d)
            rewrite_angle_line(os.path.join(d, "illumina.in"), elev, azim)
            wall = run(binary, d, [], os.path.join(d, "run.log"), args.timeout)
            t_single += wall
            print("single run %d (%g %g): %.2f s" % (k + 1, elev, azim, wall))
            single.append(d)

        # 2. one looped run
        d = os.path.join(work, "looped")
        copy_inputs(case, d)
        angles = os.path.join(d, "angles.txt")
        with open(angles, "w") as f:
            f.write("# elevation_deg azimuth_deg (geographic)\n\n")
            for elev, azim in POINTINGS:
                f.write("%g %g\n" % (elev, azim))
            f.write("\n# end\n")
        t_loop = run(binary, d, ["illumina.in", basenm + ".out", "angles.txt"], os.path.join(d, "run.log"), args.timeout)
        print("looped run (%d pointings): %.2f s" % (len(POINTINGS), t_loop))

        # 3. byte comparison of every pointing
        combined_parts = []
        for k, (elev, azim) in enumerate(POINTINGS):
            root = pointing_root(basenm, elev, azim)
            for ext in ("_result.txt", "_pcl.bin"):
                a = os.path.join(single[k], basenm + ext)
                b = os.path.join(d, root + ext)
                if not os.path.exists(b):
                    ok = False
                    print("MISSING  %s" % b)
                    continue
                same = read_bytes(a) == read_bytes(b)
                ok &= same
                print("%-8s %s == single %d %s" % ("ok" if same else "MISMATCH", root + ext, k + 1, basenm + ext))
            combined_parts.append(read_bytes(os.path.join(single[k], basenm + "_result.txt")))
            if not os.path.exists(os.path.join(d, root + ".out")):
                ok = False
                print("MISSING  %s" % os.path.join(d, root + ".out"))

        # 4. first looped pointing against reference.json
        root0 = pointing_root(basenm, *POINTINGS[0])
        result = {"summary": rr.parse_out(os.path.join(d, root0 + ".out"))}
        nbx, nby, nz = rr.read_bin(os.path.join(d, root0 + "_pcl.bin"))
        result["pcl_shape"] = [nbx, nby]
        result["pcl_nonzero"] = nz
        ref_ok, rows = rr.compare(result, ref, args.rtol, args.atol)
        ok &= ref_ok
        rr.print_table(rows, args.rtol)
        print("%-8s looped pointing 1 vs reference.json" % ("ok" if ref_ok else "MISMATCH"))

        # 5. combined record
        combined = read_bytes(os.path.join(d, basenm + "_results.txt"))
        expected = b"\n".join(combined_parts)
        same = combined == expected
        ok &= same
        print("%-8s %s_results.txt == concatenation of the single _result.txt" % ("ok" if same else "MISMATCH", basenm))

        # 6. timing
        print("timing: 3 single runs %.2f s, 1 looped run %.2f s (ratio %.2f)" % (t_single, t_loop, t_loop / t_single if t_single else float("nan")))

        # 7. nomaps: same _result.txt files, no _pcl.bin
        dn = os.path.join(work, "nomaps")
        copy_inputs(case, dn)
        shutil.copy(angles, dn)
        t_nomaps = run(binary, dn, ["illumina.in", basenm + ".out", "angles.txt", "nomaps"], os.path.join(dn, "run.log"), args.timeout)
        print("nomaps looped run (%d pointings): %.2f s" % (len(POINTINGS), t_nomaps))
        for elev, azim in POINTINGS:
            root = pointing_root(basenm, elev, azim)
            a = os.path.join(d, root + "_result.txt")
            b = os.path.join(dn, root + "_result.txt")
            same = os.path.exists(b) and read_bytes(a) == read_bytes(b)
            ok &= same
            print("%-8s nomaps %s == looped" % ("ok" if same else "MISMATCH", root + "_result.txt"))
        pcl = [n for n in os.listdir(dn) if n.endswith("_pcl.bin")]
        ok &= not pcl
        print("%-8s nomaps wrote no _pcl.bin%s" % ("ok" if not pcl else "MISMATCH", "" if not pcl else ": " + " ".join(pcl)))
        same = read_bytes(os.path.join(dn, basenm + "_results.txt")) == combined
        ok &= same
        print("%-8s nomaps %s_results.txt == looped" % ("ok" if same else "MISMATCH", basenm))
    finally:
        if args.keep:
            print("work directory kept:", work)
        else:
            shutil.rmtree(work, ignore_errors=True)

    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
