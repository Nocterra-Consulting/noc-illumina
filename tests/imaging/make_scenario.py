#!/usr/bin/env python3
"""Build the synthetic "imaging" scenario for a full-hemisphere ILLUMINA run.

Only numpy is required (no GDAL, no scipy, no illum package). The script
reproduces the pre-processing of ``illum/inputs.py`` and
``illum/inventory.py`` (from_lamps) with numpy only and writes every file
that the Fortran kernel reads:

    scenario/
      inventory.txt          invented lamp inventory (x_m y_m lumens obsth
                             obstd obstf altlp spct lop), one lamp per line
      angles.lst             1210 pointings "elev_deg azim_deg"
      summary.txt            lamp counts, lumens, reflectance, aerosol choice
      shared/                terrain / lamp height / obstacle / origin bins
      wl_<centre>/           one kernel working directory per wavelength bin:
                             illumina.in, aerosol.txt, layer.txt,
                             MolecularAbs.txt, sky_fctem_00N.dat,
                             sky_lumlp_00N.bin and symlinks into ../shared/

Domain, lamps and terrain
-------------------------
* 128 x 128 cells of 50 m (6.4 km). Row 0 is north, column 0 is west.
* Terrain: 0 m plus one Gaussian hill (250 m, sigma 600 m) centred 1.5 km
  north-east of the observer. Seen from the observer the hill top is at
  about 10 deg elevation toward azimuth 45 deg.
* Observer: domain centre, 0.2 m above ground (nalmp default).
* Type 1 = HPS (HPS_Helios.spct + 5_pcUPLIGHT.lop): a 15 x 10 "town" grid
  of 150 lamps, one lamp per cell, centred 1 km south-west of the observer.
* Type 2 = LED 4000 K (4LED_LED-4000K-Philips.spct + 1_pcUPLIGHT.lop):
  80 lamps every 50 m along an east-west road 2 km south of the observer.
* ``--lamp-scale`` thins both sets evenly. The default 0.5 keeps every
  second lamp (75 HPS on a checkerboard, 40 LED every 100 m) because the
  kernel time per pointing scales with the lamp count: with all 230
  lamps one pointing takes about 7 s (4 threads) and the 1210-pointing
  hemisphere would take about 140 min; with 115 lamps about 3.5 s and
  about 75 min. The domain size has little effect (96 cells: -12 %).
* Lamp height 8 m. Obstacles on lamp cells: height 6 m, distance 15 m,
  fill 0.5; zero elsewhere. The kernel turns obstd = 0 into
  drefle = dx (the cell size) and obsH = 0 into "no obstacle"
  (kernel/illumina.f, the obsth/obstd/obstf reads and the angmin tests).
  Upstream (inputs.py, end) would instead nearest-neighbour fill the
  obstacle arrays from the lamp cells over the whole domain.

fctem normalisation (illum/AngularPowerDistribution.py)
-------------------------------------------------------
``APD.from_txt`` reads the .lop file (columns: value, zenith angle 0..180,
0 = up). ``normalize`` divides by ``sum(vertical_profile(integrated=True))``
where the integrated profile is ``p * 2 pi * diff(-cos(mids(angles)))``,
i.e. the integral of p over the sphere with mid-point solid-angle weights
becomes 1. ``interpolate(step=1)`` is the identity for a 1-degree file.
``inventory.from_lamps`` writes ``fctem_wl_<wl>_lamp_<lop>.dat`` with two
columns (value, angle). The kernel reads the first column of each of the
181 lines and normalises again (pvalto sum, kernel/illumina.f near
"reading photometry files"), so the file normalisation only has to be
consistent. Here ``lop_norm`` does the same sum (identical to
``illum.pytools.LOP_norm``).

lumlp normalisation (illum/inputs.py + illum/inventory.py)
----------------------------------------------------------
``norm_spectrum = photopic.dat / max * 683.002`` (lm/W) on the photopic
wavelength grid ``wav`` (273..899.5 nm, 0.5 nm). Every lamp spectrum is
interpolated on ``wav`` and divided by ``trapz(S * norm_spectrum, wav)`` so
that the spectrum integrates to 1 lm. In ``from_lamps`` the spectrum of
every lamp is multiplied by its lumens and summed per cell; the lumlp
value of a bin is ``mean(spectrum[bin mask])``, the mean spectral power in
W/nm inside the bin (mask ``wav >= lo & wav < hi``). The kernel multiplies
by the bandwidth.

Reflectance
-----------
Upstream (inputs.py, "Interpolating reflectance") weights the ASTER
files by the ``reflectance`` block of the parameters (asphalt 0.8, grass
0.2) and takes the mean per bin. The same is done here, but the ASTER
curves are clamped to their edge values outside their wavelength range
(asphalt starts at 420 nm; the upstream interpolator returns 0 there).

Aerosols (illum/OPAC.py)
------------------------
Upstream builds the phase function of the aerosol mixture ("MC" here:
waso 1500, ssam 20, sscm 3.2e-3 per Aerosol_optics/combination_types,
relative humidity 80 %) from the OPAC tables. The extinction and
scattering coefficients are interpolated linearly in wavelength. Upstream
interpolates the phase function with a cubic 2-D spline; here it is
interpolated linearly between the two tabulated OPAC wavelengths that
bracket the bin centre (the "nearest available wavelength" columns of
the table) and linearly in angle onto 0..180 deg. The layer file uses the
"CC" mixture, as the nalmp default; the layer AOD is 0 so it is inert.
"""
import argparse
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
# The nalmp Lights directory (spectra, photometry, ASTER reflectance,
# photopic curve). Override with --lights or the ILLUM_LIGHTS variable.
DEFAULT_LIGHTS = os.environ.get("ILLUM_LIGHTS") or os.path.join(
    os.path.dirname(REPO),
    "artificial-light-modelling",
    "nalmp",
    "data",
    "default_folder",
    "Lights",
)

BASENM = "sky"
ANGLES = np.arange(181)
_trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy < 2.0 has only trapz

LAMP_TYPES = [
    # (name, spct file, lop file, lumens per lamp)
    ("HPS", "HPS_Helios.spct", "5_pcUPLIGHT.lop", 10000.0),
    ("LED4000K", "4LED_LED-4000K-Philips.spct", "1_pcUPLIGHT.lop", 8000.0),
]

# nalmp user_params.yml defaults
LAMBDA_MIN, LAMBDA_MAX, N_BINS = 330.0, 780.0, 5
AOD, ANGSTROM, HAER = 0.11, 0.7, 2000
PRESSURE = 101.3
DIRECT_FOV = 5.0
REFL_RADIUS = 9.99
STOP_LIMIT = 5000.0
OBS_HEIGHT = 0.2
RH = 80
REFLECTANCE_MIX = {"asphalt": 0.8, "grass": 0.2}
AEROSOL_PROFILE = "MC"
LAYER_TYPE = "CC"

OPAC_TYPES = ["inso", "waso", "soot", "ssam", "sscm", "minm", "miam",
              "micm", "mitr", "suso", "fogr"]
OPAC_RH_TYPES = ["waso", "ssam", "sscm", "suso"]
OPAC_MIX = {
    "MC": {"waso": 1500, "ssam": 20, "sscm": 3.2e-3},
    "CC": {"inso": 0.15, "waso": 2600},
}


# ----------------------------------------------------------------------
# helpers shared with tests/regression/make_case.py
# ----------------------------------------------------------------------
def save_bin(filename, data):
    """Write a 2-D array in the ILLUMINA binary format (same as illum)."""
    data = np.asarray(data, dtype=np.float32)
    ny, nx = data.shape
    head = np.array([8, nx, ny, 8], dtype=np.uint32)
    filler = np.full(data.size, 5.6e-45, dtype=np.float32)
    body = np.stack([filler, data.ravel(), filler], axis=1).ravel()
    with open(filename, "wb") as f:
        head.tofile(f)
        body.tofile(f)


def lop_norm(angles_deg, x):
    """Normalise x so that its integral over the sphere is 1.

    Same as illum.pytools.LOP_norm and APD.normalize."""
    a = np.deg2rad(angles_deg)
    mids = np.concatenate([[a[0]], 0.5 * (a[1:] + a[:-1]), [a[-1]]])
    sinx = 2 * np.pi * (np.cos(mids[:-1]) - np.cos(mids[1:]))
    return x / np.sum(x * sinx)


def input_line(values, comment, n_space=30):
    value_str = " ".join(str(v) for v in values)
    return "%-*s ! %s" % (n_space, value_str, comment)


# ----------------------------------------------------------------------
# photometry and spectra
# ----------------------------------------------------------------------
def load_lop(filename):
    """Return the 181-value normalised profile, index = zenith angle."""
    data = np.loadtxt(filename)
    values, ang = data[:, 0], data[:, 1]
    order = np.argsort(ang)
    values, ang = values[order], ang[order]
    prof = np.interp(ANGLES, ang, values, left=0.0, right=0.0)
    return lop_norm(ANGLES, prof)


def load_xy(filename, skiprows=1):
    data = np.loadtxt(filename, skiprows=skiprows)
    x, y = data[:, 0], data[:, 1]
    x, idx = np.unique(x, return_index=True)
    return x, y[idx]


def load_spectrum(filename, wav, norm_spectrum):
    """Spectrum on `wav`, normalised to 1 lm (SPD.normalize(norm))."""
    x, y = load_xy(filename)
    s = np.interp(wav, x, y, left=0.0, right=0.0)
    return s / _trapz(s * norm_spectrum, wav)


def load_aster(filename, wav):
    x, y = load_xy(filename, skiprows=0)
    return np.interp(wav, x * 1000.0, y / 100.0)  # clamped at the edges


# ----------------------------------------------------------------------
# OPAC aerosols (numpy version of illum/OPAC.py)
# ----------------------------------------------------------------------
def read_opac(filename):
    """Return (wl_um, ext, sca, pf_angles, pf[angle, wl])."""
    props = []
    with open(filename) as f:
        lines = f.readlines()
    in_props = False
    for line in lines:
        if line.startswith("#"):
            body = line[1:].split()
            if body[:2] == ["wavelength", "ext.coef"]:
                in_props = True
                continue
            if "phase" in body and "function" in body:
                in_props = False
            if in_props and len(body) == 9:
                try:
                    props.append([float(v) for v in body])
                except ValueError:
                    pass
    props = np.array(props)
    pf = np.loadtxt(filename)
    if pf.shape[1] != props.shape[0] + 1:
        sys.exit("OPAC table %s: %d wavelengths, %d phase columns"
                 % (filename, props.shape[0], pf.shape[1] - 1))
    return props[:, 0], props[:, 1], props[:, 2], pf[:, 0], pf[:, 1:]


def opac_mixture(mix, wl_nm, rh):
    """Return (ssa, phase[181] normalised to 4 pi, (wl_lo, wl_hi) used)."""
    wl = wl_nm / 1000.0
    ext_tot = sca_tot = 0.0
    pf_tot = np.zeros(181)
    bracket = None
    for typ, n in mix.items():
        rh_t = rh if typ in OPAC_RH_TYPES else 0
        fname = os.path.join(REPO, "Aerosol_optics", "OPAC_data", "%s%02d" % (typ, rh_t))
        wls, ext, sca, pf_ang, pf = read_opac(fname)
        ext_tot += n * np.interp(wl, wls, ext)
        sca_tot += n * np.interp(wl, wls, sca)
        k = np.searchsorted(wls, wl)
        k = min(max(k, 1), len(wls) - 1)
        w = (wl - wls[k - 1]) / (wls[k] - wls[k - 1])
        pf_wl = (1 - w) * pf[:, k - 1] + w * pf[:, k]
        pf_tot += n * np.interp(ANGLES, pf_ang, pf_wl)
        bracket = (wls[k - 1] * 1000.0, wls[k] * 1000.0)
    ssa = sca_tot / ext_tot
    pf_norm = 4 * np.pi * lop_norm(ANGLES, pf_tot)
    return ssa, pf_norm, bracket


def write_phase_file(filename, ssa, pf):
    with open(filename, "w") as f:
        f.write("%s # single scatering albedo\n" % ssa)
        f.write("ScatAngle PhaseFct\n")
        np.savetxt(f, np.stack([ANGLES, pf], 1), fmt="%g")


# ----------------------------------------------------------------------
# parameter file
# ----------------------------------------------------------------------
def write_illumina_in(filename, dx, x_obs, y_obs, wl, bw, refl, ntype):
    lines = [
        input_line([""], "Input file for ILLUMINA"),
        input_line([BASENM], "Root file name"),
        input_line([dx, dx], "Cell size along X [m] ; Cell size along Y [m]"),
        input_line(["aerosol.txt"], "Aerosol optical cross section file"),
        input_line(["layer.txt", 0.0, 1.0, 2000],
                   "Layer file ; Layer AOD at 500nm ; Layer angstrom ; Layer scale height [m]"),
        input_line([1], "Double scattering activated"),
        input_line([1], "Single scattering activated"),
        input_line(["%g" % wl, "%g" % bw], "Wavelength [nm] ; Bandwidth [nm]"),
        input_line(["%.6g" % refl], "Reflectance"),
        input_line([PRESSURE], "Ground level pressure [kPa]"),
        input_line([AOD, ANGSTROM, HAER],
                   "Aerosol optical depth at 500nm ; Angstrom exponent ; Aerosol scale height [m]"),
        input_line([ntype], "Number of source types"),
        input_line(["%g" % STOP_LIMIT], "Contribution threshold"),
        input_line([""], ""),
        input_line([x_obs, y_obs, OBS_HEIGHT],
                   "Observer X position ; Observer Y position ; Observer elevation above ground [m]"),
        input_line([0], "Obstacles around observer"),
        input_line([30.0, 45.0], "Elevation viewing angle ; Azimuthal viewing angle"),
        input_line([DIRECT_FOV], "Direct field of view"),
        input_line([""], ""),
        input_line([""], ""),
        input_line([""], ""),
        input_line([REFL_RADIUS], "Radius around light sources where reflextions are computed"),
        input_line([0, 0, 0], "Cloud model (0=clear) ; Cloud base altitude [m] ; Cloud fraction"),
        input_line([""], ""),
    ]
    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")


# ----------------------------------------------------------------------
# angles grid (nalmp user_params.yml defaults)
# ----------------------------------------------------------------------
def segments_to_angles(start, segs):
    out = [float(start)]
    cur = float(start)
    for inc, end in segs:
        n = int(round((end - cur) / inc))
        out.extend(np.linspace(cur, end, n + 1)[1:].tolist())
        cur = float(end)
    return np.array(out)


def default_grid():
    azi = segments_to_angles(-180, [(10, -30), (2.5, 30), (10, 180)])
    ele = segments_to_angles(0, [(2.5, 30), (5, 60), (10, 90)])
    return azi, ele


# ----------------------------------------------------------------------
# lamp inventory
# ----------------------------------------------------------------------
def build_inventory(n, dx, scale):
    """Return a list of (x_m, y_m, lumens, obsth, obstd, obstf, altlp, type).

    Coordinates are metres from the observer (x east, y north). `scale`
    thins both lamp sets (1.0 = full counts)."""
    lamps = []
    obs = (6.0, 15.0, 0.5, 8.0)
    # Town: 15 columns x 10 rows, one lamp per cell, centre 1 km SW.
    cx, cy = -707.0, -707.0
    cols = np.arange(15) - 7
    rows = np.arange(10) - 4
    for j in rows:
        for i in cols:
            lamps.append((cx + i * dx, cy - j * dx, LAMP_TYPES[0][3]) + obs + (0,))
    # Road: 80 lamps every cell on an east-west line 2 km south.
    ry = -2000.0
    for i in np.arange(80) - 40:
        lamps.append((i * dx, ry, LAMP_TYPES[1][3]) + obs + (1,))
    if scale < 1.0:
        keep = []
        for t in range(len(LAMP_TYPES)):
            idx = [k for k, l in enumerate(lamps) if l[-1] == t]
            n_keep = max(1, int(round(len(idx) * scale)))
            sel = np.linspace(0, len(idx) - 1, n_keep).round().astype(int)
            keep.extend(idx[s] for s in sel)
        lamps = [lamps[k] for k in sorted(set(keep))]
    half = n * dx / 2.0
    for l in lamps:
        if abs(l[0]) >= half or abs(l[1]) >= half:
            sys.exit("lamp at (%g, %g) m falls outside the %d-cell domain" % (l[0], l[1], n))
    return lamps


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(HERE, "scenario"))
    ap.add_argument("--size", type=int, default=128, help="domain size in cells")
    ap.add_argument("--dx", type=float, default=50.0, help="cell size [m]")
    ap.add_argument("--lamp-scale", type=float, default=0.5,
                    help="fraction of the lamps to keep; 1.0 = 150 HPS + 80 LED, "
                         "default 0.5 (75 + 40) keeps the hemisphere run under 90 min")
    ap.add_argument("--lights", default=DEFAULT_LIGHTS,
                    help="directory with the .spct/.lop/.aster/photopic.dat files "
                         "(default: $ILLUM_LIGHTS or %s)" % DEFAULT_LIGHTS)
    args = ap.parse_args(argv)

    n, dx = args.size, args.dx
    out = args.out
    shared = os.path.join(out, "shared")
    os.makedirs(shared, exist_ok=True)
    lights = args.lights
    for fname in ["photopic.dat", "asphalt.aster", "grass.aster"] + [
        t[1] for t in LAMP_TYPES] + [t[2] for t in LAMP_TYPES]:
        if not os.path.isfile(os.path.join(lights, fname)):
            sys.exit("missing %s in %s" % (fname, lights))

    c = n // 2
    # cell centre coordinates relative to the observer cell (row 0 = north)
    xs = (np.arange(n) - c) * dx
    ys = (c - np.arange(n)) * dx
    X, Y = np.meshgrid(xs, ys)

    # --- terrain -------------------------------------------------------
    hill_d = 1500.0 / np.sqrt(2.0)
    topo = 250.0 * np.exp(-((X - hill_d) ** 2 + (Y - hill_d) ** 2) / (2 * 600.0**2))
    topo = topo.astype(np.float32)

    # --- lamps ---------------------------------------------------------
    lamps = build_inventory(n, dx, args.lamp_scale)
    altlp = np.zeros((n, n), np.float32)
    obsth = np.zeros((n, n), np.float32)
    obstd = np.zeros((n, n), np.float32)
    obstf = np.zeros((n, n), np.float32)
    origin = np.zeros((n, n), np.float32)
    lumens_cell = np.zeros((len(LAMP_TYPES), n, n))
    with open(os.path.join(out, "inventory.txt"), "w") as f:
        f.write("# x_m y_m lumens obsth_m obstd_m obstf altlp_m spct lop  "
                "(x east, y north, metres from the observer)\n")
        for x, y, lum, oh, od, of, hl, t in lamps:
            f.write("%g %g %g %g %g %g %g %s %s\n" % (
                x, y, lum, oh, od, of, hl,
                LAMP_TYPES[t][1].split("_")[0], LAMP_TYPES[t][2].split("_")[0]))
            col = int(round(x / dx)) + c
            row = c - int(round(y / dx))
            if (col, row) == (c, c):
                sys.exit("a lamp falls on the observer cell")
            lumens_cell[t, row, col] += lum
            # lumen-weighted averages as in from_lamps; constant here
            altlp[row, col] = hl
            obsth[row, col] = oh
            obstd[row, col] = od
            obstf[row, col] = of

    # --- photometry ----------------------------------------------------
    fctem = [load_lop(os.path.join(lights, t[2])) for t in LAMP_TYPES]

    # --- spectra -------------------------------------------------------
    wav, phot = load_xy(os.path.join(lights, "photopic.dat"))
    norm_spectrum = phot / phot.max() * 683.002
    spct = [load_spectrum(os.path.join(lights, t[1]), wav, norm_spectrum) for t in LAMP_TYPES]

    limits = np.linspace(LAMBDA_MIN, LAMBDA_MAX, N_BINS + 1)
    bins = np.stack([limits[:-1], limits[1:]], 1)
    centres = bins.mean(1)
    bws = bins[:, 1] - bins[:, 0]
    masks = [(wav >= lo) & (wav < hi) for lo, hi in bins]

    # --- reflectance ---------------------------------------------------
    refl_curve = sum(load_aster(os.path.join(lights, k + ".aster"), wav) * w
                     for k, w in REFLECTANCE_MIX.items()) / sum(REFLECTANCE_MIX.values())
    refl = [float(np.mean(refl_curve[m])) for m in masks]

    # --- shared files --------------------------------------------------
    x_obs, y_obs = c + 1, n - c
    save_bin(os.path.join(shared, BASENM + "_topogra.bin"), topo)
    save_bin(os.path.join(shared, BASENM + "_altlp.bin"), altlp)
    save_bin(os.path.join(shared, BASENM + "_obsth.bin"), obsth)
    save_bin(os.path.join(shared, BASENM + "_obstd.bin"), obstd)
    save_bin(os.path.join(shared, BASENM + "_obstf.bin"), obstf)
    save_bin(os.path.join(shared, "origin.bin"), origin)
    shutil.copy(os.path.join(REPO, "Molecular_optics", "MolecularAbs.txt"),
                os.path.join(shared, "MolecularAbs.txt"))
    shared_files = [BASENM + s for s in ["_topogra.bin", "_altlp.bin", "_obsth.bin",
                                         "_obstd.bin", "_obstf.bin"]] + [
        "origin.bin", "MolecularAbs.txt"]

    # --- per wavelength ------------------------------------------------
    summary = []
    summary.append("Scenario: %dx%d cells, dx=%g m, observer cell (x=%d, y=%d) Fortran, %g m above ground"
                   % (n, n, dx, x_obs, y_obs, OBS_HEIGHT))
    summary.append("Hill: 250 m, sigma 600 m, centre (%.0f, %.0f) m from observer; top elevation seen from the observer ~%.1f deg"
                   % (hill_d, hill_d, np.degrees(np.arctan(topo.max() / 1500.0))))
    for t, (name, s, l, lum) in enumerate(LAMP_TYPES):
        nl = sum(1 for lp in lamps if lp[-1] == t)
        summary.append("Type %d %s: %s + %s, %d lamps x %g lm = %g lm, %d cells"
                       % (t + 1, name, s, l, nl, lum, nl * lum, int((lumens_cell[t] > 0).sum())))
    summary.append("Bins (nm): " + ", ".join("%g-%g (centre %g)" % (lo, hi, cc) for (lo, hi), cc in zip(bins, centres)))
    summary.append("Reflectance per bin (asphalt 0.8 / grass 0.2 ASTER mean): "
                   + ", ".join("%.4f" % r for r in refl))
    for b, (wl, bw) in enumerate(zip(centres, bws)):
        wdir = os.path.join(out, "wl_%g" % wl)
        os.makedirs(wdir, exist_ok=True)
        for fname in shared_files:
            link = os.path.join(wdir, fname)
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(os.path.join("..", "shared", fname), link)
        for t in range(len(LAMP_TYPES)):
            np.savetxt(os.path.join(wdir, "%s_fctem_%03d.dat" % (BASENM, t + 1)),
                       np.stack([fctem[t], ANGLES], 1))
            lumlp = lumens_cell[t] * np.mean(spct[t][masks[b]])  # W/nm
            save_bin(os.path.join(wdir, "%s_lumlp_%03d.bin" % (BASENM, t + 1)), lumlp)
        ssa, pf, br = opac_mixture(OPAC_MIX[AEROSOL_PROFILE], wl, RH)
        write_phase_file(os.path.join(wdir, "aerosol.txt"), ssa, pf)
        ssa_l, pf_l, br_l = opac_mixture(OPAC_MIX[LAYER_TYPE], wl, RH)
        write_phase_file(os.path.join(wdir, "layer.txt"), ssa_l, pf_l)
        write_illumina_in(os.path.join(wdir, "illumina.in"), dx, x_obs, y_obs,
                          wl, bw, refl[b], len(LAMP_TYPES))
        tot = [float(lumens_cell[t].sum() * np.mean(spct[t][masks[b]])) for t in range(len(LAMP_TYPES))]
        summary.append("wl_%g: aerosol %s RH%d from OPAC columns %g/%g nm (ssa %.4f); layer %s; "
                       "reflectance %.4f; lumlp totals W/nm: %s"
                       % (wl, AEROSOL_PROFILE, RH, br[0], br[1], ssa, LAYER_TYPE, refl[b],
                          ", ".join("type %d %.4g" % (t + 1, v) for t, v in enumerate(tot))))

    # --- angles --------------------------------------------------------
    azi, ele = default_grid()
    azi_geo = np.where(azi < 0, azi + 360.0, azi)
    with open(os.path.join(out, "angles.lst"), "w") as f:
        f.write("# elevation_deg azimuth_deg (geographic, 0 = north, 90 = east)\n")
        f.write("# nalmp default grid: %d azimuths x %d elevations\n" % (len(azi), len(ele)))
        for e in ele:
            for a in azi_geo:
                f.write("%g %g\n" % (e, a))
    summary.append("Angles: %d azimuths x %d elevations = %d pointings (azimuth -180 and 180 both map to 180)"
                   % (len(azi), len(ele), len(azi) * len(ele)))

    with open(os.path.join(out, "summary.txt"), "w") as f:
        f.write("\n".join(summary) + "\n")
    print("\n".join(summary))
    print("Scenario written to", out)


if __name__ == "__main__":
    main()
