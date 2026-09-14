#!/usr/bin/env bash
# Run the imaging scenario over the full hemisphere with bin/illumina_omp.
#
#   tests/imaging/run_hemisphere.sh [--probe-only] [--force] [--threads N]
#                                   [--probe "ELEV AZIM"] [--limit-min M]
#
# 1. Probe: one pointing per wavelength directory (all five in parallel,
#    OMP_NUM_THREADS threads each) with OUTFILE "probe.out". The projected
#    wall time is (slowest probe) x (number of pointings in angles.lst).
#    If it exceeds --limit-min (default 90) the script stops with exit
#    code 2 unless --force is given. Regenerate a smaller scenario with
#       python3 tests/imaging/make_scenario.py --size 96
#    or  python3 tests/imaging/make_scenario.py --lamp-scale 0.5
#    and run again.
# 2. Full run: bin/illumina_omp illumina.in sky.out ../angles.lst in every
#    wl_* directory in parallel, log in wl_*/run.log, combined results in
#    wl_*/sky_results.txt (the output root is OUTFILE without ".out").
#    The wall time is printed and stored in
#    scenario/timing.txt.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
SCEN="$HERE/scenario"
BIN="$REPO/bin/illumina_omp"
THREADS=4
PROBE_PT="20 45"
LIMIT_MIN=90
FORCE=0
PROBE_ONLY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --force) FORCE=1 ;;
    --probe-only) PROBE_ONLY=1 ;;
    --threads) THREADS=$2; shift ;;
    --probe) PROBE_PT=$2; shift ;;
    --limit-min) LIMIT_MIN=$2; shift ;;
    *) echo "unknown option $1" >&2; exit 1 ;;
  esac
  shift
done

[ -x "$BIN" ] || { echo "missing $BIN (run: make openmp)" >&2; exit 1; }
[ -f "$SCEN/angles.lst" ] || { echo "missing $SCEN (run: python3 tests/imaging/make_scenario.py)" >&2; exit 1; }
WLDIRS=$(ls -d "$SCEN"/wl_* | sort -t_ -k2 -n)
NPTS=$(grep -cve '^\s*#' -e '^\s*$' "$SCEN/angles.lst")
export OMP_NUM_THREADS=$THREADS

echo "binary: $BIN"
echo "threads per run: $THREADS, wavelength dirs in parallel: $(echo "$WLDIRS" | wc -l), pointings: $NPTS"

# ---------------------------------------------------------------- probe
echo "== probe: one pointing ($PROBE_PT) per wavelength"
echo "$PROBE_PT" > "$SCEN/probe_angles.lst"
for d in $WLDIRS; do
  ( cd "$d" && rm -f probe* && s=$(date +%s.%N) && "$BIN" illumina.in probe.out ../probe_angles.lst > probe.log 2>&1;
    e=$(date +%s.%N); printf '%s\n' "$(echo "$e - $s" | bc)" > probe.time ) &
done
wait
SLOWEST=0
for d in $WLDIRS; do
  t=$(cat "$d/probe.time")
  ok=$([ -f "$d/probe_result.txt" ] && grep -q diffuse_radiance "$d/probe_result.txt" && echo ok || echo FAILED)
  printf '  %-8s %8.2f s  %s\n' "$(basename "$d")" "$t" "$ok"
  if [ "$ok" != ok ]; then echo "probe failed in $d, see probe.log" >&2; tail -5 "$d/probe.log" >&2; exit 1; fi
  SLOWEST=$(echo "if ($t > $SLOWEST) $t else $SLOWEST" | bc)
done
PROJ_MIN=$(echo "$SLOWEST * $NPTS / 60" | bc -l)
printf 'projected full run: %.2f s/pointing x %d pointings = %.1f min (limit %d min)\n' "$SLOWEST" "$NPTS" "$PROJ_MIN" "$LIMIT_MIN"
if [ "$(echo "$PROJ_MIN > $LIMIT_MIN" | bc)" = 1 ] && [ $FORCE = 0 ]; then
  echo "projected time exceeds the limit: regenerate a smaller scenario (make_scenario.py --size 96 or --lamp-scale 0.5) or use --force" >&2
  exit 2
fi
[ $PROBE_ONLY = 1 ] && exit 0

# ------------------------------------------------------------- full run
echo "== full run: $NPTS pointings per wavelength"
START=$(date +%s)
for d in $WLDIRS; do
  ( cd "$d" && rm -f sky_e*_a* sky_results.txt && "$BIN" illumina.in sky.out ../angles.lst > run.log 2>&1; echo $? > run.exit ) &
done
wait
END=$(date +%s)
WALL=$((END - START))
{
  echo "full run wall time: $WALL s ($(echo "$WALL / 60" | bc -l | cut -c1-6) min), $NPTS pointings x $(echo "$WLDIRS" | wc -l) wavelengths"
  echo "per pointing (wall / pointings, 5 wavelengths in parallel): $(echo "$WALL / $NPTS" | bc -l | cut -c1-6) s"
  for d in $WLDIRS; do
    n=$(grep -c '^diffuse_radiance=' "$d/sky_results.txt" 2>/dev/null || echo 0)
    err=$(grep -ci 'error' "$d/run.log")
    printf '  %-8s exit %s  results %s/%s  error lines %s\n' "$(basename "$d")" "$(cat "$d/run.exit")" "$n" "$NPTS" "$err"
  done
} | tee "$SCEN/timing.txt"
