#!/usr/bin/env python3
"""Run an ILLUMINA kernel over a list of pointings and collect the results.

The script works with two kinds of kernel:

- ``--per-pointing``: a kernel without command line arguments and without
  the angle loop (the ``master`` branch). Every pointing gets its own work
  directory with symbolic links to the input files and a copy of
  ``illumina.in`` whose line 17 holds the pointing. One process per
  pointing, ``--jobs`` processes in parallel.
- default: a kernel with the angle loop (``ai-update``). The angle list is
  split into ``--jobs`` chunks; every chunk runs as one process
  ``illumina illumina.in run.out <chunk.lst>`` in its own work
  directory. ``--threads`` sets ``OMP_NUM_THREADS`` for every process.

Both modes parse the six summary numbers from every ``.out`` file (3
significant digits, the only precision the ``master`` kernel writes) and
count the line of sight steps (``Progression along the line of sight``).
The results are written to ``<out>/results.txt`` as ``key=value`` blocks
in the format of ``<basenm>_results.txt``, plus ``los_steps=``. The
wall time of the run goes to ``<out>/timing.txt``.

Example (master, one worktree build):

    python3 tests/imaging/run_pointings.py --per-pointing \\
        --binary /tmp/illum-master/bin/illumina \\
        --case tests/imaging/scenario/wl_555 \\
        --angles tests/imaging/scenario/angles_200.lst \\
        --out tests/imaging/scenario/cmp_master --jobs 20
"""
import argparse
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "tests", "regression"))
from run_regression import LABELS, parse_out  # noqa: E402

KEYS = [
    "direct_irradiance_sources",
    "direct_irradiance_reflection",
    "direct_radiance_sources",
    "direct_radiance_reflection",
    "cloud_radiance",
    "diffuse_radiance",
]

# Input files the kernel reads from the case directory (no outputs).
INPUT_SUFFIXES = ("_topogra.bin", "_altlp.bin", "_obsth.bin", "_obstd.bin",
                  "_obstf.bin")
INPUT_FILES = ("MolecularAbs.txt", "aerosol.txt", "layer.txt", "origin.bin")


def read_angles(path):
    angles = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            e, a = ln.split()[:2]
            angles.append((float(e), float(a)))
    return angles


def angle_tag(v):
    """Kernel file name tag of an angle: one decimal, '.'->'p', '-'->'m'."""
    return ("%.1f" % v).replace("-", "m").replace(".", "p")


def read_basenm_wl(case):
    with open(os.path.join(case, "illumina.in")) as f:
        lines = f.readlines()
    return lines[1].split()[0], float(lines[7].split()[0])


def link_inputs(case, basenm, work):
    os.makedirs(work, exist_ok=True)
    names = list(INPUT_FILES) + [basenm + s for s in INPUT_SUFFIXES]
    names += [n for n in os.listdir(case)
              if n.startswith(basenm + "_lumlp_") or n.startswith(basenm + "_fctem_")]
    for n in names:
        src = os.path.join(case, n)
        dst = os.path.join(work, n)
        if os.path.isfile(src) and not os.path.exists(dst):
            os.symlink(os.path.abspath(src), dst)


def write_parfile(case, work, pointing=None):
    with open(os.path.join(case, "illumina.in")) as f:
        lines = f.readlines()
    if pointing is not None:
        lines[16] = "%-30s ! Elevation viewing angle ; Azimuthal viewing angle\n" % (
            "%g %g" % pointing)
    with open(os.path.join(work, "illumina.in"), "w") as f:
        f.writelines(lines)


def los_steps(path):
    with open(path, errors="replace") as f:
        return sum(1 for ln in f if "Progression along the line of sight" in ln)


def run(cmd, cwd, env, log):
    with open(log, "w") as f:
        proc = subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT, env=env)
    if proc.returncode != 0:
        raise RuntimeError("kernel failed (%d) in %s, see %s" % (proc.returncode, cwd, log))


def collect(out_path, elev, azim, wl):
    rec = {"elevation_deg": elev, "azimuth_deg": azim, "wavelength_nm": wl}
    vals = parse_out(out_path)
    for lab, key in zip(LABELS, KEYS):
        rec[key] = vals[lab]
    rec["los_steps"] = los_steps(out_path)
    return rec


def write_results(path, records):
    with open(path, "w") as f:
        for r in records:
            for k in ["elevation_deg", "azimuth_deg", "wavelength_nm"] + KEYS:
                f.write("%s=%14.6E\n" % (k, r[k]))
            f.write("los_steps=%d\n\n" % r["los_steps"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--case", required=True, help="kernel directory with illumina.in and the inputs")
    ap.add_argument("--angles", required=True, help="angles list (elevation azimuth per line)")
    ap.add_argument("--out", required=True, help="work/output directory")
    ap.add_argument("--jobs", type=int, default=20, help="parallel processes (default 20)")
    ap.add_argument("--threads", type=int, default=1, help="OMP_NUM_THREADS per process (default 1)")
    ap.add_argument("--per-pointing", action="store_true",
                    help="one process per pointing with line 17 replaced (kernel without angle loop)")
    ap.add_argument("--keep-logs", action="store_true", help="keep the kernel stdout logs")
    args = ap.parse_args(argv)

    binary = os.path.abspath(args.binary)
    case = os.path.abspath(args.case)
    out = os.path.abspath(args.out)
    angles = read_angles(args.angles)
    basenm, wl = read_basenm_wl(case)
    os.makedirs(out, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads))

    t0 = time.time()
    records = [None] * len(angles)
    if args.per_pointing:
        def one(k):
            elev, azim = angles[k]
            work = os.path.join(out, "p_e%s_a%s" % (angle_tag(elev), angle_tag(azim)))
            link_inputs(case, basenm, work)
            write_parfile(case, work, (elev, azim))
            run([binary], work, env, os.path.join(work, "run.log"))
            records[k] = collect(os.path.join(work, basenm + ".out"), elev, azim, wl)
            if not args.keep_logs:
                os.remove(os.path.join(work, "run.log"))

        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            list(ex.map(one, range(len(angles))))
    else:
        nchunk = min(args.jobs, len(angles))
        chunks = [list(range(i, len(angles), nchunk)) for i in range(nchunk)]

        def one(c):
            work = os.path.join(out, "chunk_%02d" % c)
            link_inputs(case, basenm, work)
            write_parfile(case, work)
            lst = os.path.join(work, "angles.lst")
            with open(lst, "w") as f:
                for k in chunks[c]:
                    f.write("%g %g\n" % angles[k])
            run([binary, "illumina.in", "run.out", lst], work, env, os.path.join(work, "run.log"))
            for k in chunks[c]:
                elev, azim = angles[k]
                if len(chunks[c]) == 1:
                    name = "run.out"
                else:
                    name = "run_e%s_a%s.out" % (angle_tag(elev), angle_tag(azim))
                records[k] = collect(os.path.join(work, name), elev, azim, wl)
            if not args.keep_logs:
                os.remove(os.path.join(work, "run.log"))

        with ThreadPoolExecutor(max_workers=nchunk) as ex:
            list(ex.map(one, range(nchunk)))
    wall = time.time() - t0

    write_results(os.path.join(out, "results.txt"), records)
    mode = "per-pointing" if args.per_pointing else "angle-loop"
    line = ("%s  binary=%s  mode=%s  pointings=%d  jobs=%d  threads=%d  wall=%.1f s  (%.2f s per pointing per process)\n"
            % (time.strftime("%Y-%m-%d %H:%M:%S"), binary, mode, len(angles), args.jobs, args.threads,
               wall, wall * min(args.jobs, len(angles)) / len(angles)))
    with open(os.path.join(out, "timing.txt"), "a") as f:
        f.write(line)
    print(line.strip())
    print("results:", os.path.join(out, "results.txt"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
