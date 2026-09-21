# Synthetic regression harness

This directory holds four small synthetic ILLUMINA cases:

- `case_small`: flat terrain, clear sky.
- `case_hill`: one Gaussian hill, clear sky (see "Hill case").
- `case_cloud_hill`: the hill terrain under a cloud (see "Cloudy
  cases").
- `case_cloud_2nd`: flat terrain under a high cloud with a steep view
  (see "Cloudy cases").

Each case runs in a few seconds and does not need the `illum` Python
package or GDAL. Only `python3` with `numpy` is required.

The check compares:

- the six summary numbers at the end of `<basenm>.out` (printed by the
  kernel with 7 significant digits, format `E14.7E2`),
- the non-zero pixels of the contribution map `<basenm>_pcl.bin`
  (full float32 precision).

## Files

- `make_case.py` writes every input file that the kernel reads into a
  case directory (default `case_small/`).
- `run_regression.py` runs a kernel binary in the case directory and
  compares the results with `case_small/reference.json`.
- `run_angles.py` checks the multi-pointing loop of the kernel (see
  "Multi-pointing check" below).
- `reference.json` in each of the four case directories is the stored
  reference (committed). All other files in `case_*/` are generated and
  ignored by git.

## Kernel command line

The kernel accepts four optional positional arguments:

    illumina [PARFILE [OUTFILE [ANGLES_FILE [maps|nomaps]]]]

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
- `maps` or `nomaps`: `maps` (the default when the argument is absent)
  writes the per-pointing contribution map `<root>_pcl.bin`; `nomaps`
  skips that file (3 MB per pointing) and writes every other output
  (`.out`, `_result.txt`, `_results.txt`). Any other value stops the
  run with a message. The regression harness runs with the default, so
  the references hold the maps; `run_angles.py` checks that `nomaps`
  gives the same `_result.txt` files and no `_pcl.bin`.

All arguments are copied as is, so a path with `/`, a space or a
comma works. The kernel prints the resolved file names and the number
of pointings at start-up. A file name that the kernel builds is at most
512 characters (`maxnam`); a longer name stops the run with a message
that gives the variable and the length needed. The base name (line 2 of
the parameter file) and every input name derived from it (`_topogra`,
`_obsth`, `_obstd`, `_altlp`, `_obstf`, `_fctem_NNN`, `_lumlp_NNN`) share
that limit; the former 72-character limit on these names is gone (#73).

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

## Near-field exclusion radius

The kernel drops every source-voxel pair closer than an exclusion
radius (10 m by default, the historical value hard-wired in `omemax`).
An optional 25th line of the parameter file sets that radius in metres;
a file without the line keeps 10 m and gives bit-identical results.
When the cell size is at most twice the exclusion radius the kernel
prints a `WARNING` line to stdout and to the `.out` header: most of the
near-field signal is then discarded, and cells smaller than the
exclusion radius are outside the valid range of the model. A case built
with `make_case.py --dx 5` triggers the warning; the four regression
cases (100 m cells) do not.

## Run the check

    make test

or

    make
    python3 tests/regression/run_regression.py --binary bin/illumina
    python3 tests/regression/run_regression.py --binary bin/illumina \
        --case tests/regression/case_hill
    python3 tests/regression/run_regression.py --binary bin/illumina \
        --case tests/regression/case_cloud_hill
    python3 tests/regression/run_regression.py --binary bin/illumina \
        --case tests/regression/case_cloud_2nd
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
6. reports the wall time of the three single runs and of the looped run,
7. runs the looped run once more with the fourth argument `nomaps` and
   asserts that every `_result.txt` and the combined `synth_results.txt`
   are byte-identical to the looped run and that no `_pcl.bin` exists.

The script accepts `--binary`, `--case`, `--rtol`, `--atol` and
`--timeout` like `run_regression.py`. Use `--keep` to keep the
temporary work directories. The exit status is 0 on PASS and 1 on
FAIL.

## OpenMP check

    make test-omp

or

    make openmp
    python3 tests/regression/run_omp.py --serial bin/illumina --omp bin/illumina_omp

`make test` runs this check on all four cases after the serial checks
when `bin/illumina_omp` exists. The kernel has one `!$omp parallel do`
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
- `case_hill` exposed the observer horizon bug (#59). Before the fix
  the serial kernel read `zhoriz` in the line of sight test
  (`angze1-zhoriz.lt.0.00001`) after the source and reflection loops
  had overwritten it with the horizon of the last processed cell. The
  OpenMP build keeps `zhoriz` private in the loop, so the 1-thread run
  saw the observer horizon and differed from the serial run (14 line
  of sight steps serial, 13 with OpenMP). The kernel now stores the
  observer horizon in `zhorob`; serial and 1-thread runs are
  byte-identical on both cases. `case_small` (flat terrain) never
  showed the difference.

## Regenerate the case

    rm -rf tests/regression/case_small
    python3 tests/regression/make_case.py

`make_case.py` accepts `--size N` (default 64), `--dx M` (default 100 m),
`--pad 512` (embed the domain in a 512 x 512 zero array as the Python
side does), `--double-scattering 0|1`, `--stop-limit X`, `--lamps N`
(about N extra lamps on a regular grid, for timing cases; default 0)
and `--cloud MODEL BASE_M FRACTION` (line 23 of `illumina.in`:
`cloudt`, `cloudbase`, `cloudfrac`; default `0 0 0`, clear sky).
`MODEL` is 0 for a clear sky, 1 thin cirrus, 2 thick cirrus, 3
altostratus/altocumulus, 4 cumulus/cumulonimbus, 5 stratocumulus.

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

`case_hill/` is a second case with terrain. Its `reference.json` is
committed and `make test` runs it after `case_small`. The harness
regenerates the inputs from `make_case_args` when they are missing; to
generate them by hand:

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

## Comparing against `main`

`main` is the pre-refactor kernel (`master` merged with the DUG
`updated_batches` build). Three fixes separate it from `ai-update`, and
between them they account for every difference the four cases showed
before the reflection disc change (#90, see "Reference history"), which
now also moves the reflected channels on every case. If a
comparison against `main` turns up anything outside this table, it is a
finding and not an expected change.

| case | expected difference against `main` | cause | detail |
|------|------------------------------------|-------|--------|
| `case_small` | none; the contribution map is bit-identical and the six numbers agree to the 3 significant figures `main` prints | — | — |
| `case_hill` | clear sky, but the contribution map moves: the pixel behind the hill falls from 3.7196695e-03 to 1.0552055e-06, and the three bright pixels rise by 0.37 % | observer horizon (#59) | "Hill case" and "Reference history" below |
| `case_cloud_hill` | diffuse radiance -78.8 %; the map moves between x0.028 and x14.9 | cloud double count (upstream `e63f3e0`) | `case_cloud_hill` below |
| `case_cloud_2nd` | cloud radiance -88.6 % and diffuse -95.9 % against `main` (the second fix accounts for -82.4 % of it, on top of the first) | both cloud fixes | `case_cloud_2nd` below |

The cleanest observable of terrain behaviour is the line-of-sight step
count, not any radiance. The harness records it as `los_steps` (with the
furthest horizontal distance as `los_max_dist_m`) and compares it
exactly, so a change in where the sight line stops is reported before
any number moves. On a real scenario the same signal was decisive: at
0.4 and 0.8 degrees elevation `main` ran the sight line out to 2.2 and
2.5 km, into the sources, where the fixed kernel stops at 30 m.

Two traps when comparing:

- `main` prints the summary with `E10.3E2`, three significant figures,
  where `ai-update` prints `E14.7E2`. Comparing raw values invents
  relative differences of up to 0.5 % on quantities that are provably
  bit-identical. Round the `ai-update` side to three figures, or compare
  `_pcl.bin`, which carries full precision on both.
- The contribution map moves much further than the summary numbers
  wherever a fix applies. A summary-level comparison understates the
  change for anything downstream that reads the map.

The observer horizon fix is the one to watch on real terrain. Its error
is not signed: `main` can report too much flux past a ridge, as it does
here, or too little, depending on which source cell the loop happened to
process last. It only ever touches pointings between the ground and the
local horizon, in the azimuth sector the relief occupies.

## Cloudy cases

Two cases run with a cloud layer (`cloudt` other than 0). The suite had
no cloudy case before, so the cloud code path was never checked. The two
upstream cloud fixes below changed the result of both cases.

Generate them by hand with:

    python3 tests/regression/make_case.py --out tests/regression/case_cloud_hill \
        --hill 200 500 1060.7 1060.7 --hill-lamps --obs-height 0.2 \
        --view 30 45 --cloud 3 1000 100
    python3 tests/regression/make_case.py --out tests/regression/case_cloud_2nd \
        --size 64 --view 80 45 --cloud 3 5000 100

### `case_cloud_hill` (double counting of the cloud radiance)

The hill case with an altostratus/altocumulus layer (model 3) at a base
of 1000 m and a cloud fraction of 100 %. The viewing elevation rises
from 3 deg to 30 deg, so the line of sight clears the summit and reaches
the cloud base. The domain, the terrain and the eight lamps are the same
as `case_hill`.

The case covers the "cloud base detection in the line of sight" fix
(upstream `e63f3e0`): the cloud intensity `icloud` entered the result
twice, once through `isourc=isourc+icloud` and once through `fctcld`.
The direct channels and the cloud radiance stay the same; the diffuse
radiance and the contribution map change.

| quantity | before the fix | after the fix | change |
|----------|---------------:|--------------:|-------:|
| irdirect | 1.552939e-04 | 1.552939e-04 | 0 |
| irrdirect | 9.309442e-05 | 9.309442e-05 | 0 |
| cloud | 4.082893e-05 | 4.082893e-05 | 0 |
| diffuse | 2.013058e-04 | 4.271468e-05 | -78.8 % |
| pcl(37,36) | 2.2812013e-02 | 3.3937421e-01 | x 14.9 |
| pcl(37,38) | 4.3630410e-02 | 2.5877854e-01 | x 5.9 |
| pcl(39,36) | 6.7651004e-02 | 2.8230387e-01 | x 4.2 |
| pcl(51,50) | 2.5098833e-01 | 2.5995538e-02 | x 0.10 |
| pcl(52,23) | 2.5450671e-01 | 7.1637584e-03 | x 0.028 |
| pcl(46,48) | 1.3889244e-01 | 3.8579978e-02 | x 0.28 |
| pcl(49,45) | 2.1796846e-01 | 4.1151989e-02 | x 0.19 |
| pcl(24,51) | 3.5507083e-03 | 6.6520600e-03 | x 1.9 |

The second cloud fix (upstream `1c46ee1`) leaves this case unchanged:
with a cloud base of 1000 m no second order scattering voxel reaches the
cloud. The kernel now says so: it prints one `WARNING: the second-order
scattering volume lies entirely above the cloud base` line per pointing
to stdout and to the `.out` file when every second order scattering cell
above ground is discarded by the cloud base (#77). The line appears on
this case only; the clear-sky cases and `case_cloud_2nd` (base 5000 m)
do not print it, and the references are unchanged.

### `case_cloud_2nd` (wrong flux in the cloud radiance after a first scattering)

Flat terrain, 64 x 64 cells of 100 m, the default six lamps of
`case_small`, an altostratus/altocumulus layer (model 3) at a base of
5000 m, a cloud fraction of 100 % and a viewing elevation of 80 deg.
This is the only configuration found that reaches the
reflection-then-second-scattering cloud block: the steep view keeps the
line of sight inside the domain up to the cloud base, and the high base
leaves room for a second scattering voxel below the cloud.

The case covers both fixes. The first fix removes the double counting;
the second fix replaces `fldif2` (the flux at the scattering voxel) by
`fdif2` (the flux at the line of sight voxel) in the cloud term, which
is the flux that illuminates a cloud on the sight line.

| quantity | before the fixes | after fix 1 | after fix 1 + fix 2 |
|----------|-----------------:|------------:|--------------------:|
| irdirect | 1.095854e-04 | 1.095854e-04 | 1.095854e-04 |
| irrdirect | 3.908474e-09 | 3.908474e-09 | 3.908474e-09 |
| cloud | 7.764644e-06 | 7.764644e-06 | 8.815393e-07 |
| diffuse | 3.553758e-05 | 8.349037e-06 | 1.465932e-06 |

Fix 1 alone lowers the diffuse radiance by 76.5 %. Fix 2 then lowers the
cloud radiance by 88.6 % and the diffuse radiance by a further 82.4 %.

### Serial and OpenMP

The four cloud proximity tests of the first fix relax from
`cloudbase-z_c .le. iz*scal` to `.le. 1.20*iz*scal`. Neither new case
measures that part: a build with the old test gives the same cloud
radiance and the same diffuse radiance on both cases. The other parts of
the first fix carry the whole change.

The OpenMP region keeps the cloud path on one thread (`if(cloudt.eq.0)`
clause), because `icloud` is a running sum that the loop also reads. The
two cloudy cases are therefore a check that the serial and the OpenMP
builds still agree. They do: the 1-thread run is byte-identical to the
serial run and the 4-thread run has a maximum relative difference of
0.0e+00 on both cases.

## Reference history

- All four `reference.json` were rebuilt on 2026-09-22 from the serial
  kernel after the reflection disc change (#90). Before the change the
  reflection box was `boxx=nint(reflsiz/dx)` cells wide, a step
  function of `dx` (one cell above `dx=2*reflsiz`, nine cells from
  `2/3*reflsiz` to `2*reflsiz`, 25 below), and only the source cell
  used the sub-grid side `reflsiz`; every neighbour reflected with its
  full area `dx*dy`. At `dx=10` and `reflsiz=9.99` the reflecting area
  was 99.8 + 8 x 100 = 900 m2 for a parameter that asks for a 10 m
  radius. The kernel now treats the reflecting region as the disc of
  radius `reflsiz` centred on the source: the box covers the disc
  (`boxx=ceiling(reflsiz/dx)`), and each cell in the box reflects with
  the fraction of its area inside the disc (`discfrac`, analytic), as a
  square of side `dx*sqrt(afrac)`. The source cell follows the same
  rule, so the total reflecting area is `pi*reflsiz**2` for every cell
  size. Checks made with the change:

  - `case_small` at `dx=100`, `reflsiz` swept over 5, 9.99, 30, 60, 90,
    160 (serial build, `irrdirect` / `diffuse`):

    | reflsiz | old irrdirect | new irrdirect | old diffuse | new diffuse |
    |--------:|--------------:|--------------:|------------:|------------:|
    | 5 | 5.372469e-09 | 1.687611e-08 | 8.985140e-05 | 9.017894e-05 |
    | 9.99 | 2.144338e-08 | 6.730727e-08 | 9.029627e-05 | 9.098759e-05 |
    | 30 | 1.928516e-07 | 6.023938e-07 | 9.167710e-05 | 9.219757e-05 |
    | 60 | 1.435750e-05 | 2.319441e-06 | 9.245966e-05 | 9.250508e-05 |
    | 90 | 1.529045e-05 | 4.886324e-06 | 9.261908e-05 | 9.251531e-05 |
    | 160 | 3.422712e-05 | 1.405974e-05 | 9.276698e-05 | 9.263289e-05 |

    The old `irrdirect` jumps by a factor 74 between 30 and 60 (the
    box opens from one cell to nine) and again at 160 (25 cells). The
    new `irrdirect` is smooth and monotone; `irrdirect/reflsiz**2` goes
    6.8e-10, 6.7e-10, 6.7e-10, 6.4e-10, 6.0e-10, 5.5e-10 (a slow fall
    from the lamp photometry, no step).
  - Flat twin of the hill case (`--hill 0 500 1060.7 1060.7
    --hill-lamps --obs-height 1.5 --view 30 45`, eight lamps placed in
    metres so the lamp set is the same at every cell size) at
    `dx=100, 50, 25, 10` (`--size 6400/dx`). `irrdirect` at a fixed
    `reflsiz`, old -> new:

    | reflsiz | dx=100 | dx=50 | dx=25 | dx=10 |
    |--------:|-------:|------:|------:|------:|
    | 9.99 old | 1.023e-04 | 1.101e-04 | 1.134e-04 | 3.203e-04 |
    | 9.99 new | 2.119e-04 | 2.281e-04 | 2.349e-04 | 1.832e-04 |
    | 30 old | 3.192e-04 | 3.796e-04 | 4.263e-04 | 4.109e-04 |
    | 30 new | 4.003e-04 | 4.205e-04 | 3.760e-04 | 3.904e-04 |
    | 60 old | 4.258e-04 | 4.605e-04 | 4.413e-04 | 4.311e-04 |
    | 60 new | 4.475e-04 | 4.428e-04 | 4.364e-04 | 4.272e-04 |

    The reflected part of the diffuse radiance (`diffuse(reflsiz)`
    minus `diffuse(0.001)`), new kernel: 2.42e-06, 2.36e-06, 2.42e-06,
    1.90e-06 at 9.99; 4.57e-06, 4.36e-06, 3.88e-06, 4.04e-06 at 30;
    5.11e-06, 4.59e-06, 4.50e-06, 4.42e-06 at 60. The old kernel gave
    1.17e-06, 1.14e-06, 1.17e-06, 3.32e-06 at 9.99 (x2.8 between
    `dx=100` and `dx=10`). The new result is constant across `dx` to
    within 13 % at 9.99, 6 % at 30 and 2.5 % at 60; the remaining
    spread is the cell-centre discretisation of a 10 m high lamp. The
    contribution per unit `reflsiz**2` is not constant across
    `reflsiz` on this case: the lamps are 10 m high with a downward
    photometry, so the ground irradiance saturates beyond about 10 m.
  - Size of the change on the four cases (old reference -> new
    reference, serial build). `irdirect` and `direct` are unchanged on
    every case.

    | case | quantity | old | new | change |
    |------|----------|----:|----:|-------:|
    | `case_small` | irrdirect | 2.1443380e-08 | 6.7307270e-08 | x 3.14 |
    | | diffuse | 9.0296270e-05 | 9.0987590e-05 | +0.8 % |
    | | los_steps | 22 | 24 | stop criterion moved |
    | | pcl(29,36) | 5.5025902e-04 | 1.0481671e-03 | +90 % |
    | | pcl(31,28) | 3.8446407e-04 | 7.4386160e-04 | +93 % |
    | | pcl(34,38) | 2.1928316e-03 | 3.9066863e-03 | +78 % |
    | | pcl(36,30) | 1.2836505e-03 | 2.3565716e-03 | +84 % |
    | | pcl(37,36) | 9.9020159e-01 | 9.8269689e-01 | -0.8 % (mast) |
    | | pcl(38,35) | 5.3871297e-03 | 9.2478590e-03 | +72 % |
    | `case_hill` | irrdirect | 1.0371790e-04 | 2.1556250e-04 | x 2.08 |
    | | rdirect | 7.2969160e-03 | 1.5147860e-02 | x 2.08 |
    | | diffuse | 9.2855650e-05 | 1.1834080e-04 | +27 % |
    | | pcl(24,51) | 4.9881542e-06 | 7.9156462e-06 | +59 % |
    | | pcl(37,36) | 9.1420996e-01 | 9.1044128e-01 | -0.4 % |
    | | pcl(37,38) | 4.1023631e-02 | 4.2822793e-02 | +4.4 % |
    | | pcl(39,36) | 4.4753052e-02 | 4.6715781e-02 | +4.4 % |
    | | pcl(46,48) | 9.1573844e-07 | 1.1874182e-06 | +30 % |
    | | pcl(49,45) | 9.7678719e-07 | 1.2665809e-06 | +30 % |
    | | pcl(51,50) | 1.0552055e-06 | 1.3565627e-06 | +29 % |
    | | pcl(52,23) | 5.3718609e-06 | 8.5245429e-06 | +59 % |
    | `case_cloud_hill` | irrdirect | 9.3094420e-05 | 1.9348400e-04 | x 2.08 |
    | | cloud | 4.0828930e-05 | 6.5228840e-05 | +60 % |
    | | diffuse | 4.2714680e-05 | 6.8425990e-05 | +60 % |
    | | pcl(24,51) | 6.6520600e-03 | 7.4333455e-03 | +12 % |
    | | pcl(37,36) | 3.3937421e-01 | 3.3871546e-01 | -0.2 % |
    | | pcl(37,38) | 2.5877854e-01 | 2.5817788e-01 | -0.2 % |
    | | pcl(39,36) | 2.8230387e-01 | 2.8164867e-01 | -0.2 % |
    | | pcl(46,48) | 3.8579978e-02 | 3.8508479e-02 | -0.2 % |
    | | pcl(49,45) | 4.1151989e-02 | 4.1075714e-02 | -0.2 % |
    | | pcl(51,50) | 2.5995538e-02 | 2.6435262e-02 | +1.7 % |
    | | pcl(52,23) | 7.1637584e-03 | 8.0051431e-03 | +12 % |
    | `case_cloud_2nd` | irrdirect | 3.9084740e-09 | 1.2268060e-08 | x 3.14 |
    | | cloud | 8.8153930e-07 | 1.3183160e-06 | +50 % |
    | | diffuse | 1.4659320e-06 | 2.1926830e-06 | +50 % |
    | | pcl(29,36) | 1.1128492e-01 | 1.2377754e-01 | +11 % |
    | | pcl(31,28) | 1.0164575e-01 | 1.1356333e-01 | +12 % |
    | | pcl(34,38) | 1.4267236e-01 | 1.5757380e-01 | +10 % |
    | | pcl(36,30) | 2.1995874e-01 | 2.4195480e-01 | +10 % |
    | | pcl(37,36) | 2.4120682e-01 | 1.6128626e-01 | -33 % (mast) |
    | | pcl(38,35) | 1.8323134e-01 | 2.0184425e-01 | +10 % |

    Why the numbers move: every case uses `reflsiz=9.99` at `dx=100`,
    so the old kernel reflected from one 9.99 m square (99.8 m2) and
    the new kernel reflects from a disc of 313.5 m2 (x 3.14). On the
    flat cases with a far, high mast lamp the reflected irradiance
    scales with that area (x 3.14). On the hill cases the lamps are
    10 m high and the reflecting square sits right under them, so the
    solid angle grows less than the area (x 2.08). The ground lamps
    gain weight in the contribution map and the mast lamp loses it.
    `los_steps` on `case_small` moves from 22 to 24 because the
    `1/stoplim` stop criterion compares each step with the flux
    accumulated so far, which is now larger. The 1-thread OpenMP run
    stays byte-identical to the serial run on every case and the
    4-thread run agrees within 1.2e-07.
- The `.out` summary format changed from `E10.3E2` to `E14.7E2`
  (upstream `e63f3e0`), so `case_small/reference.json` and
  `case_hill/reference.json` were rebuilt on 2026-09-18 from the patched
  serial kernel. The refresh is a precision change only:
  - every one of the twelve stored numbers, rounded back to the old
    `E10.3E2` form, gives exactly the old stored value (for example
    `case_hill` diffuse 9.2855650e-05 -> 0.929E-04 = 9.29e-05, direct
    1.0575880e-02 -> 0.106E-01 = 1.06e-02),
  - `pcl_shape` and every entry of `pcl_nonzero` are unchanged, and
    `synth_pcl.bin` is byte-identical to the pre-patch run for both
    clear-sky cases.
- `case_cloud_hill/reference.json` and `case_cloud_2nd/reference.json`
  (2026-09-18) are new. They are built from the patched serial kernel
  (both cloud fixes applied); there is no earlier reference to preserve.
- `case_small/reference.json` is built from `main` (e33b0f4, `master` 56745d2 merged with `updated_batches`). It is bit-identical to the earlier `master` reference.
- `case_hill/reference.json` (2026-09-14) is built from the `ai-update`
  serial kernel with the observer horizon fix (#59, variable `zhorob`).
  The `git_commit` field records the `ai-update` checkout at run time
  (3416aa6 before the fix was committed). A reference built from `main`
  (e33b0f4) in a worktree was bit-identical to the pre-fix `ai-update`
  serial kernel; the flat `case_small` reference is unchanged by the
  fix. Before/after numbers of the fix on `case_hill` (serial build,
  `synth.out` and `synth_pcl.bin`):

  | quantity | main / pre-fix | fixed | note |
  |----------|---------------:|------:|------|
  | line of sight steps | 14 | 13 | last step was inside the terrain |
  | irdirect | 1.720e-04 | 1.720e-04 | unchanged |
  | irrdirect | 1.040e-04 | 1.040e-04 | unchanged |
  | direct | 1.060e-02 | 1.060e-02 | unchanged |
  | rdirect | 7.300e-03 | 7.300e-03 | unchanged |
  | diffuse | 9.320e-05 | 9.290e-05 | -0.3 % |
  | pcl(37,36) | 9.1080999e-01 | 9.1420996e-01 | +0.4 % (in front of the hill) |
  | pcl(37,38) | 4.0871121e-02 | 4.1023631e-02 | +0.4 % |
  | pcl(39,36) | 4.4586673e-02 | 4.4753052e-02 | +0.4 % |
  | pcl(24,51) | 5.0313579e-06 | 4.9881542e-06 | -0.9 % (flank) |
  | pcl(52,23) | 5.4183879e-06 | 5.3718609e-06 | -0.9 % (flank) |
  | pcl(46,48) | 1.0199171e-06 | 9.1573844e-07 | -10 % (behind the hill) |
  | pcl(49,45) | 1.0879110e-06 | 9.7678719e-07 | -10 % (behind the hill) |
  | pcl(51,50) | 3.7196695e-03 | 1.0552055e-06 | -100 % (behind the hill; the stale step lit it) |

  The fixed serial run is byte-identical to the pre-fix 1-thread OpenMP
  run (`synth_result.txt` and `synth_pcl.bin`), and to the fixed
  1-thread OpenMP run. The flat twin of the case (`--hill 0 500 1060.7
  1060.7 --hill-lamps --obs-height 0.2 --view 3 45`) computes 30 line
  of sight steps, diffuse 3.03e-05 and behind-the-hill pcl weights of
  0.06 to 0.17, so the hill changes the result materially.
