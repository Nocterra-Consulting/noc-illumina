# Synthetic regression harness

This directory holds a small synthetic ILLUMINA case. The case runs in
a few seconds and does not need the `illum` Python package or GDAL. Only
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

    illumina [PARFILE [OUTFILE [ANGLES_FILE]]]

- `PARFILE`: the parameter file (default `illumina.in`).
- `OUTFILE`: the output file name (default `<basenm>.out`, where
  `basenm` is line 2 of the parameter file). Input file names always
  use `basenm`. The output root `<root>` is `OUTFILE` without a
  trailing `.out`.
- `ANGLES_FILE`: a list of pointings, one `elevation_deg azimuth_deg`
  pair per line (geographic azimuth, the convention of line 17 of the
  parameter file). Blank lines and lines that start with `#` are
  skipped. Every elevation must be within -90..90 deg. The kernel stops
  with a message when the file is missing, has no pointing, or has a
  line that cannot be parsed. Without this argument the single pointing
  of the parameter file is used.

Arguments 1 and 2 are read with a list-directed read, so a path that
contains `/` is cut at the slash. Argument 3 is copied as is, so any
path works. The kernel prints the resolved file names and the number of
pointings at start-up.

The kernel loads the domain once and loops over the pointings. Output
files:

- one pointing: `<root>.out`, `<root>_result.txt` and
  `<basenm>_pcl.bin` (unchanged names; the contribution map keeps the
  `basenm` root as before),
- several pointings: `<root>_e<elev>_a<azim>.out`, `_pcl.bin` and
  `_result.txt`, where the angles are written with one decimal, `.`
  replaced by `p` and a leading `-` replaced by `m` (for example
  `synth_e30p0_a45p0.out`, `synth_em5p0_a350p0_result.txt`,
  `synth_e0p0_a0p0_pcl.bin`). Two pointings that round to the same
  tenth of a degree get the same file names.
- always: `<root>_results.txt` with the `key=value` block of every
  pointing, blocks separated by a blank line. The file is overwritten
  on every run.

## Run the check

    make test

or

    make
    python3 tests/regression/run_regression.py --binary bin/illumina
    python3 tests/regression/run_angles.py --binary bin/illumina

The harness regenerates the inputs when `<case>/illumina.in` is
missing, with the `make_case_args` stored in `reference.json`. Options
of `run_regression.py`:

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
   (`30 45` = the reference pointing, `10 120`, `60 300`),
2. runs the binary once with an angles file that lists the three
   pointings (`illumina illumina.in synth.out angles.txt`),
3. asserts that every `_result.txt` and every `_pcl.bin` of the looped
   run is byte-identical to the matching single run,
4. asserts that the first pointing of the looped run matches
   `reference.json`,
5. asserts that the combined `synth_results.txt` equals the three
   single `_result.txt` files joined with one blank line,
6. reports the wall time of the three single runs and of the looped run.

The script accepts `--binary`, `--case`, `--rtol`, `--atol` and
`--timeout` like `run_regression.py`. Use `--keep` to keep the
temporary work directories. The exit status is 0 on PASS and 1 on
FAIL.

## OpenMP check

    make test-omp

or

    make openmp
    python3 tests/regression/run_omp.py --serial bin/illumina --omp bin/illumina_omp

`make test` runs this check after the two serial checks when
`bin/illumina_omp` exists. The kernel has one `!$omp parallel do`
region: the loop over the source cells (`x_s`, `y_s`) inside each line
of sight step and each lamp type in the scattered-radiance section of
`kernel/illumina.f`. The serial build (`bin/illumina`) treats the
directives as comments and is unchanged.

`run_omp.py` copies the inputs of a case into temporary work
directories and:

1. runs the serial binary,
2. runs the OpenMP binary with `OMP_NUM_THREADS=1` and asserts that
   `<basenm>_result.txt` and `<basenm>_pcl.bin` are byte-identical to
   the serial run,
3. runs the OpenMP binary with every other thread count of `--threads`
   (default `1,4`) and asserts that the six results of
   `<basenm>_result.txt` and the non-zero pixels of `<basenm>_pcl.bin`
   agree with the serial run within `--rtol` (default 1e-4; the
   reduction changes the summation order, so bit-identity is not
   expected),
4. prints the wall time and speed-up of every run and the largest
   relative difference of every threaded run.

Options: `--serial`, `--omp`, `--case`, `--threads 1,2,4`, `--rtol`,
`--atol`, `--timeout`, `--keep` (keep the work directories) and
`--verbose` (print every compared value). The exit status is 0 on PASS
and 1 on FAIL.

For a timing run on the large case:

    python3 tests/regression/run_omp.py --case tests/regression/case_large \
        --threads 1,2,4

Notes:

- With clouds on (`cloudt` other than 0) the region runs on one thread
  (`if(cloudt.eq.0)` clause): `icloud` is a running sum that the loop
  also reads, so the result depends on the order of the source cells.
- OMPFLAGS in the Makefile adds `-fmax-stack-var-size=65536`.
  `-fopenmp` implies `-frecursive`, which puts every local array on the
  stack; `zondif(3000000,3)` (36 MB) then overflows the default 8 MB
  stack and the kernel stops with a segmentation fault. The flag keeps
  arrays above 64 KB in static storage (main program only) and every
  smaller array on the stack. gfortran prints a warning about this
  override for every source file; the warning is expected.
- The check on `case_hill` is known to differ at 1 thread. The serial
  code reads `zhoriz` in the line of sight test
  (`angze1-zhoriz.lt.0.00001`) after the source loop has overwritten
  it with the horizon of the last source or reflecting cell. The
  OpenMP build keeps `zhoriz` private in the loop, so the test sees the
  observer horizon as intended. The serial run computes 14 line of
  sight steps, the 1-thread OpenMP run 13. A serial build that stores
  the observer horizon in its own variable is byte-identical to the
  OpenMP run. `case_small` (flat terrain) is unaffected.

## Regenerate the case

    rm -rf tests/regression/case_small
    python3 tests/regression/make_case.py

`make_case.py` accepts `--size N` (default 64), `--dx M` (default 100 m),
`--pad 512` (embed the domain in a 512 x 512 zero array as the Python
side does), `--double-scattering 0|1`, `--stop-limit X` and `--lamps N`
(about N extra lamps on a regular grid, for timing cases; default 0).

For a timing case with many lamps:

    python3 tests/regression/make_case.py --size 256 --lamps 400 \
        --out tests/regression/case_large

`case_large/` is ignored by git (no `reference.json`).

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

`case_hill/` is a second case with terrain. It has no committed
`reference.json` and `make test` does not run it. Generate it with:

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
