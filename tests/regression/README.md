# Synthetic regression harness

This directory holds a small synthetic ILLUMINA case. The case runs in
about 2 s and does not need the `illum` Python package or GDAL. Only
`python3` with `numpy` is required.

The check compares:

- the six summary numbers at the end of `<basenm>.out` (printed by the
  kernel with 3 significant digits, format `E10.3E2`),
- the non-zero pixels of the contribution map `<basenm>_pcl.bin`
  (full float32 precision).

## Files

- `make_case.py` writes every input file that the kernel reads into a
  case directory (default `case_small/`).
- `run_regression.py` runs a kernel binary in the case directory and
  compares the results with `case_small/reference.json`.
- `run_angles.py` checks the multi-pointing loop of the kernel (see
  "Multi-pointing check" below).
- `case_small/reference.json` is the stored reference (committed). All
  other files in `case_small/` are generated and ignored by git.

## Kernel command line

The kernel accepts three optional positional arguments:

    illumina [PARFILE [OUTROOT [ANGLES_FILE]]]

- `PARFILE`: the parameter file (default `illumina.in`).
- `OUTROOT`: the root of the output file names (default: `basenm`, line
  2 of the parameter file). Input file names always use `basenm`.
- `ANGLES_FILE`: a list of pointings, one `elevation_deg azimuth_deg`
  pair per line (geographic azimuth, the convention of line 17 of the
  parameter file). Blank lines and lines that start with `#` are
  skipped. Without this file the single pointing of the parameter file
  is used.

The kernel loads the domain once and loops over the pointings. Output
files:

- one pointing: `<OUTROOT>.out`, `<OUTROOT>_pcl.bin`,
  `<OUTROOT>_result.txt` (unchanged names),
- several pointings: `<OUTROOT>_e<elev>_a<azim>.out`, `_pcl.bin` and
  `_result.txt`, where the angles are written with one decimal, `.`
  replaced by `p` and a leading `-` replaced by `m` (for example
  `synth_e30p0_a45p0.out`, `synth_em5p0_a350p0_result.txt`),
- always: `<OUTROOT>_results.txt` with the `key=value` block of every
  pointing, blocks separated by a blank line.

## Run the check

    make test

or

    make
    python3 tests/regression/run_regression.py --binary bin/illumina

The harness regenerates the inputs when `<case>/illumina.in` is
missing, with the `make_case_args` stored in `reference.json`. Options:

    --binary PATH    kernel binary (default bin/illumina)
    --case DIR       case directory (default tests/regression/case_small)
    --rtol X         relative tolerance (default 1e-5)
    --atol X         absolute tolerance (default 0)
    --update         rewrite reference.json from this run

The exit status is 0 on PASS and 1 on FAIL.

## Multi-pointing check

    python3 tests/regression/run_angles.py --binary bin/illumina

`make test` runs this check after `run_regression.py`. The script:

1. runs the binary three times with no arguments, each time with a
   different pointing on line 17 of a temporary copy of `illumina.in`
   (`30 45`, `10 120`, `60 300`),
2. runs the binary once with an angles file that lists the three
   pointings,
3. asserts that every `_result.txt` and every `_pcl.bin` of the looped
   run is byte-identical to the matching single run,
4. asserts that the first pointing of the looped run matches
   `reference.json`,
5. asserts that the combined `_results.txt` equals the three individual
   `_result.txt` files,
6. reports the wall time of the three single runs and of the looped run.

Use `--keep` to keep the temporary work directories.

## OpenMP check

    make test-omp

or

    make openmp
    python3 tests/regression/run_omp.py

`make test` runs this check after the two serial checks when
`bin/illumina_omp` exists. The script runs the serial binary once and the
OpenMP binary once per thread count (default `--threads 1,4`). It asserts:

1. `OMP_NUM_THREADS=1`: `_result.txt` and `_pcl.bin` are byte-identical
   to the serial run,
2. more threads: the six results and the non-zero `_pcl.bin` pixels agree
   with the serial run within `--rtol` (default 1e-4). The OpenMP
   reductions change the order of the floating-point sums, so the last
   bits can differ.

The script prints the wall time of every run. Use `--keep` to keep the
output files (`omp_serial*`, `omp_t<N>*`).

The kernel parallelises the loop over the source cells inside every line
of sight step (`kernel/illumina.f`, the `x_s`/`y_s` loop of the scattered
light section). The cloud case (`cloudt.ne.0`) runs serial because
`icloud` is a running sum that the loop also reads.

For a timing case with many lamps:

    python3 tests/regression/make_case.py --size 256 --lamps 400 \
        --out tests/regression/case_large
    python3 tests/regression/run_omp.py --case tests/regression/case_large \
        --threads 1,2,4

`case_large/` is ignored by git (no `reference.json`).

## Regenerate the case

    rm -rf tests/regression/case_small
    python3 tests/regression/make_case.py

`make_case.py` accepts `--size N` (default 64), `--dx M` (default 100 m),
`--pad 512` (embed the domain in a 512 x 512 zero array as the Python
side does), `--double-scattering 0|1`, `--stop-limit X` and `--lamps N`
(about N extra lamps on a regular grid, for timing cases; default 0).

## Rebuild the reference

Build the reference from a kernel that is known to be correct, for
example the `main` branch in a separate worktree:

    git worktree add /tmp/illum-master main
    cp Makefile /tmp/illum-master/
    make -C /tmp/illum-master
    python3 tests/regression/run_regression.py \
        --binary /tmp/illum-master/bin/illumina --update
    git worktree remove /tmp/illum-master

`reference.json` records the binary path, the git commit, the date and
the wall time of the reference run.

## Hill case

`case_hill/` is a second case with terrain. `make test` runs it after
`case_small`. Regenerate it with the arguments stored in
`case_hill/reference.json` (`make_case_args`):

    python3 tests/regression/make_case.py --out tests/regression/case_hill \
        --hill 200 500 1060.7 1060.7 --hill-lamps --obs-height 0.2 --view 3 45

- Domain: 64 x 64 pixels, 100 m pixels.
- Terrain: 0 m plus one Gaussian hill, 200 m high, sigma 500 m, centred
  1.5 km north-east of the observer (1060.7 m east, 1060.7 m north). The
  terrain at the observer cell is 2.2 m; the summit is at 7.5 deg
  elevation seen from the observer.
- Observer: pixel (33, 32), 0.2 m above ground.
- Lamps: eight lamps of one type, 10 m high, placed relative to the hill
  centre (`--hill-lamps`): three in front of the hill (cells (37,36),
  (37,38), (39,36)), two far out on the flanks ((24,51) and (52,23)) and
  three behind the hill ((46,48), (49,45), (51,50)).
- Viewing: elevation 3 deg, azimuth 45 deg (toward the hill, below the
  summit), direct FOV 5 deg. The line of sight meets the terrain about
  500 m from the observer: the kernel computes 13 line of sight steps
  instead of 30 on flat terrain, and the lamps behind the hill get a
  `_pcl.bin` weight of about 1e-6 instead of 0.06 to 0.17.
- Atmosphere: as `case_small`.

`make_case.py` options for this case: `--hill H SIGMA_M OFFSET_X_M
OFFSET_Y_M` (Gaussian hill added to the terrain, offsets in metres east
and north of the observer cell), `--hill-lamps` (the eight-lamp layout
above), `--obs-height Z` (default 10 m) and `--view ELEV AZIM` (default
`30 45`).

### Reference history

- `case_small/reference.json` is built from `main` (e33b0f4, `master` 56745d2 merged with `updated_batches`). It is bit-identical to the earlier `master` reference.
