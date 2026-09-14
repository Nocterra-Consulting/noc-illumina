# Imaging scenario: full-hemisphere synthetic run

This directory holds a synthetic but realistic ILLUMINA scenario that
runs the kernel over the full default sky grid of `nalmp`
(55 azimuths x 22 elevations = 1210 pointings) in 5 wavelength bins.
Only `python3` with `numpy` is required to build the inputs; the kernel
is `bin/illumina_omp` (`make openmp`).

## Files

- `make_scenario.py` writes `scenario/` (ignored by git).
- `run_hemisphere.sh` probes one pointing per wavelength, projects the
  total time, then runs the hemisphere in every `scenario/wl_*` directory
  in parallel. It calls `bin/illumina_omp illumina.in sky.out
  ../angles.lst` (kernel arguments `PARFILE OUTFILE ANGLES_FILE`; the
  output root is `OUTFILE` without `.out`, so the per-pointing files are
  `sky_e<elev>_a<azim>.out` and the combined file is `sky_results.txt`).
- `run_pointings.py` runs one kernel (with or without the angle loop)
  over an angles list and writes `results.txt`; `compare_master.py`
  compares two such result sets (see "Comparison with `master`").
- `scenario/summary.txt` lists lamp counts, lumens, reflectance per bin
  and the aerosol tables used (written by `make_scenario.py`).
- `scenario/timing.txt` and `scenario/run_hemisphere.log` hold the probe
  and full-run timings (written by `run_hemisphere.sh`).
- `scenario/wl_<centre>/sky_results.txt` are the deliverables: one
  `key=value` block per pointing (`elevation_deg`, `azimuth_deg`,
  `wavelength_nm`, `direct_irradiance_sources`,
  `direct_irradiance_reflection`, `direct_radiance_sources`,
  `direct_radiance_reflection`, `cloud_radiance`, `diffuse_radiance`).

## Scenario

- Domain: 128 x 128 cells of 50 m (6.4 km). Row 0 is north.
- Terrain: 0 m plus one Gaussian hill, 250 m high, sigma 600 m, centred
  1.5 km north-east of the observer (azimuth 45 deg). From the observer
  the hill top is at about 10 deg elevation, so the pointings at
  elevation 0 to 10 deg toward azimuth 45 deg are cut by the horizon.
- Observer: domain centre (Fortran cell x=65, y=64), 0.2 m above ground.
- Lamp type 1: HPS, `HPS_Helios.spct` + `5_pcUPLIGHT.lop`, 10 000 lm per
  lamp, a 15 x 10 "town" grid of cells centred 1 km south-west of the
  observer.
- Lamp type 2: LED 4000 K, `4LED_LED-4000K-Philips.spct` +
  `1_pcUPLIGHT.lop`, 8 000 lm per lamp, along an east-west road 2 km
  south of the observer, 4 km long.
- `--lamp-scale` thins both sets evenly. The default 0.5 keeps every
  second lamp: 75 HPS (checkerboard, 750 klm) and 40 LED (every 100 m,
  320 klm). `--lamp-scale 1` gives the full 150 HPS (1.5 Mlm) and 80 LED
  (640 klm). See "Timing" for the reason.
- Lamp height 8 m. Obstacles on lamp cells: height 6 m, distance 15 m,
  fill factor 0.5; zero on all other cells. The kernel treats
  `obstd = 0` as `drefle = dx` (cell size) and `obsH = 0` as no
  obstacle. Upstream `inputs.py` would fill the obstacle arrays from the
  lamp cells by nearest neighbour over the whole domain.
- `scenario/inventory.txt` is the invented lamp inventory:
  `x_m y_m lumens obsth obstd obstf altlp spct lop`, coordinates in
  metres from the observer (x east, y north).
- Atmosphere (nalmp `user_params.yml` defaults): clear sky, single and
  double scattering on, AOD 0.11 at 500 nm, Angstrom 0.7, aerosol scale
  height 2000 m, pressure 101.3 kPa, relative humidity 80 %, direct FOV
  5 deg, reflection radius 9.99 m, stop limit 5000, no particle layer
  (layer AOD 0), no cloud.
- Wavelength bins: 330-780 nm in 5 bins of 90 nm, centres 375, 465, 555,
  645 and 735 nm. One kernel directory `scenario/wl_<centre>/` per bin.
- Angles: `scenario/angles.lst`, the nalmp default grid. Azimuths
  -180..180 (10 deg steps to -30, 2.5 deg steps to 30, 10 deg steps to
  180) are converted to geographic 0..360 by adding 360 to negative
  values. Both -180 and 180 map to 180, so the grid contains that
  pointing twice (the same output file names, two identical blocks in
  `sky_results.txt`). Elevations 0..90 (2.5 deg steps to 30, 5 to 60,
  10 to 90).

## How the inputs reproduce the upstream pre-processing

The `illum` package cannot be imported here, so `make_scenario.py`
re-implements the relevant parts with numpy.

### fctem (angular photometry)

`illum/AngularPowerDistribution.py`: `from_txt` reads the `.lop` file
(columns value, zenith angle 0..180 with 0 = up). `normalize` divides by
`sum(vertical_profile(integrated=True))`, where the integrated profile is
`p * 2 pi * diff(-cos(mids(angles)))`; the integral of `p` over the
sphere becomes 1. `interpolate(step=1)` is the identity for a 1 deg
file. `illum/inventory.py` (`from_lamps`) writes
`fctem_wl_<wl>_lamp_<lop>.dat` with two columns (value, angle 0..180).
The kernel (`kernel/illumina.f`, "reading photometry files") reads the
first value of each of the 181 lines and normalises again with
`pvalto = sum(pval * 2 pi sin(theta) dtheta)`. `make_scenario.py` uses
the same mid-point solid-angle sum (`lop_norm`, identical to
`illum.pytools.LOP_norm`) and writes `sky_fctem_001.dat` (HPS) and
`sky_fctem_002.dat` (LED) in every `wl_*` directory.

### lumlp (lamp flux per bin)

`illum/inputs.py`: `norm_spectrum = photopic.dat / max * 683.002` lm/W on
the photopic grid `wav` (273..899.5 nm, 0.5 nm). Every `.spct` is
interpolated on `wav` and divided by `trapz(S * norm_spectrum, wav)`
(`SPD.normalize(norm)`), so the spectrum integrates to 1 lm.
`inventory.from_lamps`: per cell, `fctem = sum(spectrum * lumens)` and
the lumlp value of a bin is `mean(fctem[wav >= lo & wav < hi])`, the
mean spectral power in W/nm inside the bin; the kernel multiplies by the
bandwidth. `make_scenario.py` writes `sky_lumlp_001.bin` and
`sky_lumlp_002.bin` per bin with exactly this value.

### Reflectance

`illum/inputs.py` ("Interpolating reflectance") weights the ASTER curves
by the `reflectance` block of the parameters (asphalt 0.8, grass 0.2)
and takes the mean per bin. `make_scenario.py` does the same but clamps
the ASTER curves to their edge values outside their range (asphalt
starts at 420 nm; the upstream interpolator returns 0 below that). The
values are 0.0415, 0.0612, 0.0929, 0.0926 and 0.1660 for the five bins.

### Aerosols

`illum/OPAC.py` builds the aerosol mixture from the OPAC tables in
`Aerosol_optics/OPAC_data`. The mixture is `MC` (waso 1500, ssam 20,
sscm 3.2e-3; `Aerosol_optics/combination_types`) with the RH-80 tables
`waso80`, `ssam80`, `sscm80`. The extinction and scattering coefficients
are interpolated linearly in wavelength (as upstream). The phase function
is interpolated linearly between the two tabulated OPAC wavelengths that
bracket the bin centre (upstream uses a cubic 2-D spline) and linearly
in angle onto 0..180 deg, then normalised to 4 pi. Tables used:

| bin centre | OPAC wavelength columns | single scattering albedo |
|-----------:|------------------------:|-------------------------:|
| 375 nm | 350 / 400 nm | 0.9968 |
| 465 nm | 450 / 500 nm | 0.9975 |
| 555 nm | 550 / 600 nm | 0.9975 |
| 645 nm | 600 / 650 nm | 0.9976 |
| 735 nm | 700 / 750 nm | 0.9976 |

`layer.txt` uses the `CC` mixture (nalmp default); the layer AOD is 0 so
it has no effect.

## Regenerate and run

    make openmp
    python3 tests/imaging/make_scenario.py          # default: --lamp-scale 0.5
    tests/imaging/run_hemisphere.sh                 # probe, then full run

`run_hemisphere.sh` options: `--probe-only`, `--force` (run even if the
projection exceeds the limit), `--threads N` (default 4),
`--probe "ELEV AZIM"` (default `20 45`), `--limit-min M` (default 90).
The script runs the five `wl_*` directories in parallel with
`OMP_NUM_THREADS=4` (20 threads on a 20-core machine), logs to
`scenario/wl_*/run.log` and writes `scenario/timing.txt`.

`make_scenario.py` options: `--size N` (default 128), `--dx M` (default
50), `--lamp-scale F` (default 0.5), `--out DIR`, `--lights DIR` (the
nalmp `Lights` directory with the `.spct`, `.lop`, `.aster` and
`photopic.dat` files; default: the `ILLUM_LIGHTS` environment variable,
else `../artificial-light-modelling/nalmp/data/default_folder/Lights`
next to this repository).

Quick end-to-end check (one wavelength, three pointings, about 7 s with
4 threads):

    printf '20 45\n2.5 40\n60 180\n' > tests/imaging/scenario/angles3.lst
    ( cd tests/imaging/scenario/wl_555 && \
      OMP_NUM_THREADS=4 ../../../../bin/illumina_omp illumina.in sky.out ../angles3.lst )
    grep -c '^diffuse_radiance=' tests/imaging/scenario/wl_555/sky_results.txt   # 3

## Comparison with `master` (issue #61)

`run_pointings.py` runs a kernel over an angles list and collects the six
summary numbers of every `.out` file (3 significant digits, the only
precision the `master` kernel writes) plus the number of line of sight
steps. `--per-pointing` serves a kernel without the angle loop (one
process per pointing, line 17 of a temporary `illumina.in` replaced,
inputs linked); the default mode splits the list into `--jobs` chunks
and runs the angle loop of the current kernel. `compare_master.py`
compares two result sets pointing by pointing (any `key=value` block
file: `sky_results.txt`, `<root>_results.txt` or `results.txt` of
`run_pointings.py`).

    awk 'NR<=2 || $1==0 || $1==2.5 || $1==5 || $1==10' \
        tests/imaging/scenario/angles.lst > tests/imaging/scenario/angles_200.lst
    git worktree add /tmp/illum-master main && cp Makefile /tmp/illum-master/
    make -C /tmp/illum-master
    python3 tests/imaging/run_pointings.py --per-pointing \
        --binary /tmp/illum-master/bin/illumina --case tests/imaging/scenario/wl_555 \
        --angles tests/imaging/scenario/angles_200.lst \
        --out tests/imaging/scenario/cmp_master --jobs 20
    python3 tests/imaging/run_pointings.py --binary bin/illumina \
        --case tests/imaging/scenario/wl_555 --angles tests/imaging/scenario/angles_200.lst \
        --out tests/imaging/scenario/cmp_aiupdate --jobs 20
    python3 tests/imaging/compare_master.py tests/imaging/scenario/cmp_master/results.txt \
        tests/imaging/scenario/cmp_aiupdate/results.txt --label-a master --label-b ai-update
    git worktree remove /tmp/illum-master

`angles_200.lst` holds the 55 azimuths at elevations 0, 2.5, 5 and 10 deg
(220 lines, 216 distinct pointings because azimuth 180 occurs twice).
Wavelength bin 555 nm only.

### Results (master 56745d2 vs ai-update with the #59 fix, serial builds)

These numbers were measured on the `ai-update-rebased` branch, from
which these scripts were ported. The kernel of this branch carries the
same #59 fix (`zhorob`), so the same behaviour is expected; the full
216-pointing comparison has not been repeated here.

| quantity | compared | both zero | differing | max rel diff | median rel diff |
|----------|---------:|----------:|----------:|-------------:|----------------:|
| direct_irradiance_sources | 216 | 96 | 0 | 0 | 0 |
| direct_irradiance_reflection | 216 | 96 | 0 | 0 | 0 |
| direct_radiance_sources | 216 | 205 | 0 | 0 | 0 |
| direct_radiance_reflection | 216 | 205 | 0 | 0 | 0 |
| diffuse_radiance | 216 | 32 | 20 | 5.87e-01 | 0 |

The four direct terms are identical on every pointing. The number of
line of sight steps differs on 39 pointings, all at elevation 2.5 deg
(28) and 5 deg (11); elevations 0 and 10 deg are identical. The diffuse
radiance differs on 20 of them (relative difference
`|b - a| / max(|a|, |b|)`):

| pointing | master | ai-update | rel diff | steps master | steps ai-update |
|----------|---:|---:|---:|---:|---:|
| el 2.5 az 357.5 | 4.410e-10 | 1.820e-10 | 5.87e-01 | 21 | 10 |
| el 2.5 az 0.0 | 3.450e-10 | 1.870e-10 | 4.58e-01 | 19 | 10 |
| el 2.5 az 100.0 | 3.990e-10 | 6.210e-10 | 3.57e-01 | 18 | 27 |
| el 2.5 az 110.0 | 4.200e-10 | 6.250e-10 | 3.28e-01 | 19 | 27 |
| el 5.0 az 15.0 | 5.660e-10 | 4.010e-10 | 2.92e-01 | 26 | 16 |
| el 2.5 az 120.0 | 4.980e-10 | 6.750e-10 | 2.62e-01 | 21 | 28 |
| el 5.0 az 17.5 | 4.970e-10 | 3.760e-10 | 2.43e-01 | 25 | 15 |
| el 2.5 az 2.5 | 2.520e-10 | 1.920e-10 | 2.38e-01 | 15 | 10 |
| el 5.0 az 70.0 | 4.410e-10 | 3.380e-10 | 2.34e-01 | 25 | 14 |
| el 2.5 az 130.0 | 6.860e-10 | 8.070e-10 | 1.50e-01 | 24 | 29 |
| el 2.5 az 70.0 | 4.360e-11 | 3.800e-11 | 1.28e-01 | 21 | 3 |
| el 2.5 az 355.0 | 5.330e-10 | 6.080e-10 | 1.23e-01 | 23 | 27 |

Three groups of pointings:

1. Toward the hill (elevation 2.5 and 5 deg, azimuth 5 to 60 deg): the
   line of sight meets the terrain. `master` computes one step more than
   `ai-update` (for example 4 instead of 3 steps at el 2.5 az 40); that
   step lies inside the terrain and adds less than the last printed
   digit. 26 pointings, diffuse radiance identical at 3 digits.
2. Grazing the flanks of the hill (el 2.5 az 357.5 to 2.5 and az 70 to
   90; el 5 az 15, 17.5 and 70): `master` continues far past the true
   horizon (up to 21 steps instead of 3) and overestimates the diffuse
   radiance by 8 to 59 %.
3. Beside the hill, line of sight clear (el 2.5 az 100 to 130 and az
   347.5 to 355): `master` cuts the line of sight early (18 to 26 steps
   instead of 27 to 29) and underestimates the diffuse radiance by 3 to
   36 %.

Cause: before the #59 fix the observer horizon test at every line of
sight step compared the viewing zenith angle with `zhoriz`, which the
source and reflection loops of the previous step had overwritten with
the horizon of the last lamp cell. The stale value said "blocked" or
"clear" depending on that lamp's own horizon, not the observer's. The
pre-fix `ai-update` kernel (with the `horizon.f` walk limit changed from
512 cells to `nbx`, commit a65b46b) reproduces `master` exactly on all
216 pointings, so the walk-limit change has no effect on this scenario
and every difference above comes from the #59 fix. The OpenMP build
(`bin/illumina_omp`, 4 threads) agrees with the fixed serial build on
every pointing at 3 digits.

Timing on 20 cores (WSL2), 220 pointings of `wl_555`: `master`
per-pointing, 20 processes, 182 s; `ai-update` serial angle loop in 20
chunks, 132 s; `ai-update` OpenMP 4 threads in 5 chunks, 148 s (the
machine is saturated either way). The results and work directories
(`scenario/cmp_*`, about 50 MB each) are ignored by git.

## Timing

Measured on 20 cores (WSL2), 5 wavelengths in parallel, 4 threads each:

- Full lamp set (230 lamps): 7.0 s per pointing (elevation 20 deg),
  5.9 s at elevation 2.5 deg, 7.9 s at elevation 60 deg, projected
  142 min for 1210 pointings. A 96 x 96 domain saves only 12 % (6.2 s),
  because the kernel time scales with the number of lamp cells, not
  with the domain size.
- Half lamp set (115 lamps, the default): 3.5-3.9 s per pointing,
  projected 72-78 min. The shipped results use this set.

The measured wall time of the full run is in `scenario/timing.txt`.
