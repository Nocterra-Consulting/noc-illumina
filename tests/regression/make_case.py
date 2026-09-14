#!/usr/bin/env python3
"""Generate a small synthetic ILLUMINA case for regression testing.

Only numpy is required (no GDAL, no illum package). The generated
directory contains every file that the Fortran kernel reads:

    illumina.in                 parameter file, read with list-directed reads
    <basenm>_topogra.bin        ground elevation [m]              (2-D binary)
    <basenm>_altlp.bin          lamp height above ground [m]      (2-D binary)
    <basenm>_obsth.bin          sub-grid obstacle height [m]      (2-D binary)
    <basenm>_obstd.bin          sub-grid obstacle distance [m]    (2-D binary)
    <basenm>_obstf.bin          sub-grid obstacle fill factor 0-1 (2-D binary)
    <basenm>_lumlp_001.bin      lamp flux of source type 1        (2-D binary)
    <basenm>_fctem_001.dat      181 angular photometry values, one per line
    origin.bin                  VIIRS flag 0/1 per pixel          (2-D binary)
    aerosol.txt                 SSA, header, 181 lines "angle phase"
    layer.txt                   same format as aerosol.txt
    MolecularAbs.txt            copied from Molecular_optics/

The 2-D binary format is the gfortran unformatted sequential layout
written by illum.pytools.save_bin: a header record (nbx, nby) and then
one 4-byte record per value, row nby first, column 1 first.
"""
import argparse
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))


def save_bin(filename, data):
    """Write a 2-D array in the ILLUMINA binary format (same as illum)."""
    data = np.asarray(data, dtype=np.float32)
    ny, nx = data.shape
    # Record 1: two int32 -> length marker 8, nx, ny, length marker 8.
    head = np.array([8, nx, ny, 8], dtype=np.uint32)
    # Records 2..: one float32 each -> marker 4, value, marker 4.
    # 5.6e-45 is the float32 whose bit pattern is the integer 4.
    filler = np.full(data.size, 5.6e-45, dtype=np.float32)
    body = np.stack([filler, data.ravel(), filler], axis=1).ravel()
    with open(filename, "wb") as f:
        head.tofile(f)
        body.tofile(f)


def henyey_greenstein(g, angles_deg):
    mu = np.cos(np.deg2rad(angles_deg))
    return (1.0 - g**2) / (1.0 + g**2 - 2.0 * g * mu) ** 1.5


def lop_norm(angles_deg, x):
    """Normalise x so that its integral over the sphere is 1 (as illum)."""
    a = np.deg2rad(angles_deg)
    mids = np.concatenate([[a[0]], 0.5 * (a[1:] + a[:-1]), [a[-1]]])
    sinx = 2 * np.pi * (np.cos(mids[:-1]) - np.cos(mids[1:]))
    return x / np.sum(x * sinx)


def write_phase_file(filename, ssa, g):
    angles = np.arange(181)
    pf = 4 * np.pi * lop_norm(angles, henyey_greenstein(g, angles))
    with open(filename, "w") as f:
        f.write("%g # single scatering albedo\n" % ssa)
        f.write("ScatAngle PhaseFct\n")
        np.savetxt(f, np.stack([angles, pf], 1), fmt="%g")


def write_photometry(filename):
    """Street-light like photometry: mostly downward, 5 % uplight."""
    theta = np.deg2rad(np.arange(181))  # 0 = zenith, 180 = nadir
    down = np.clip(-np.cos(theta), 0, None)  # non-zero for theta > 90
    up = 0.05 * np.clip(np.cos(theta), 0, None)
    p = down + up
    np.savetxt(filename, p, fmt="%.6f")


def input_line(values, comment, n_space=30):
    value_str = " ".join(str(v) for v in values)
    return "%-*s ! %s" % (n_space, value_str, comment)


def write_illumina_in(filename, basenm, dx, x_obs, y_obs, args):
    lines = [
        input_line([""], "Input file for ILLUMINA"),
        input_line([basenm], "Root file name"),
        input_line([dx, dx], "Cell size along X [m] ; Cell size along Y [m]"),
        input_line(["aerosol.txt"], "Aerosol optical cross section file"),
        input_line(
            ["layer.txt", 0.0, 1.0, 2000],
            "Layer file ; Layer AOD at 500nm ; Layer angstrom ; Layer scale height [m]",
        ),
        input_line([int(args.double_scattering)], "Double scattering activated"),
        input_line([1], "Single scattering activated"),
        input_line([550.0, 10.0], "Wavelength [nm] ; Bandwidth [nm]"),
        input_line([0.2], "Reflectance"),
        input_line([101.3], "Ground level pressure [kPa]"),
        input_line(
            [0.11, 0.7, 2000],
            "Aerosol optical depth at 500nm ; Angstrom exponent ; Aerosol scale height [m]",
        ),
        input_line([1], "Number of source types"),
        input_line([args.stop_limit], "Contribution threshold"),
        input_line([""], ""),
        input_line(
            [x_obs, y_obs, args.obs_height],
            "Observer X position ; Observer Y position ; Observer elevation above ground [m]",
        ),
        input_line([0], "Obstacles around observer"),
        input_line(list(args.view), "Elevation viewing angle ; Azimuthal viewing angle"),
        input_line([5.0], "Direct field of view"),
        input_line([""], ""),
        input_line([""], ""),
        input_line([""], ""),
        input_line([9.99], "Radius around light sources where reflextions are computed"),
        input_line([0, 0, 0], "Cloud model (0=clear) ; Cloud base altitude [m] ; Cloud fraction"),
        input_line([""], ""),
    ]
    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--out",
        default=os.path.join(HERE, "case_small"),
        help="output case directory (default: tests/regression/case_small)",
    )
    ap.add_argument("--size", type=int, default=64, help="domain size in pixels")
    ap.add_argument("--dx", type=float, default=100.0, help="pixel size [m]")
    ap.add_argument(
        "--pad",
        type=int,
        default=0,
        help="pad the arrays to this width (e.g. 512) around the domain; 0 = no padding",
    )
    ap.add_argument(
        "--double-scattering",
        type=int,
        default=1,
        choices=[0, 1],
        help="enable second order scattering (default 1)",
    )
    ap.add_argument("--stop-limit", type=float, default=5000.0)
    ap.add_argument(
        "--lamps",
        type=int,
        default=0,
        help="add about N extra lamps of the same type on a regular grid (timing cases)",
    )
    ap.add_argument(
        "--hill",
        type=float,
        nargs=4,
        metavar=("H", "SIGMA_M", "OFFSET_X_M", "OFFSET_Y_M"),
        default=None,
        help="add a Gaussian hill of height H [m] and width SIGMA_M [m] centred "
        "OFFSET_X_M east and OFFSET_Y_M north of the observer cell",
    )
    ap.add_argument(
        "--hill-lamps",
        action="store_true",
        help="replace the default lamp set by 8 lamps around the hill centre "
        "(3 in front of the hill, 2 on the flanks, 3 behind it); needs --hill",
    )
    ap.add_argument("--obs-height", type=float, default=10.0, help="observer height above ground [m]")
    ap.add_argument(
        "--view",
        type=float,
        nargs=2,
        metavar=("ELEV_DEG", "AZIM_DEG"),
        default=[30.0, 45.0],
        help="viewing elevation and geographic azimuth [deg] (line 17 of illumina.in)",
    )
    ap.add_argument("--basenm", default="synth")
    args = ap.parse_args(argv)
    if args.hill_lamps and args.hill is None:
        sys.exit("--hill-lamps needs --hill")

    n = args.size
    os.makedirs(args.out, exist_ok=True)

    # Arrays are indexed [row(y), column(x)] with row 0 at the top (north).
    topo = np.zeros((n, n), dtype=np.float32)
    altlp = np.zeros((n, n), dtype=np.float32)
    obsth = np.zeros((n, n), dtype=np.float32)
    obstd = np.zeros((n, n), dtype=np.float32)
    obstf = np.zeros((n, n), dtype=np.float32)
    lumlp = np.zeros((n, n), dtype=np.float32)
    origin = np.zeros((n, n), dtype=np.float32)

    c = n // 2
    # Cell centre coordinates in metres from the observer cell (x east,
    # y north; row 0 is north).
    xm = (np.arange(n) - c) * args.dx
    ym = (c - np.arange(n)) * args.dx
    X, Y = np.meshgrid(xm, ym)

    def cell(mx, my):
        """(column, row) of the cell that holds the point (mx, my) [m]."""
        return c + int(round(mx / args.dx)), c - int(round(my / args.dx))

    if args.hill is not None:
        h, sigma, ox, oy = args.hill
        topo += (h * np.exp(-((X - ox) ** 2 + (Y - oy) ** 2) / (2.0 * sigma**2))).astype(np.float32)

    if args.hill_lamps:
        # Eight lamps of one type, 10 m high, placed relative to the hill
        # centre: u points from the observer to the hill centre, v is
        # perpendicular (to the left). Three lamps in front of the hill
        # (seen from the observer), two far out on the flanks (their own
        # horizon toward the line of sight is clear) and three behind it.
        h, sigma, ox, oy = args.hill
        d = np.hypot(ox, oy)
        u = np.array([ox, oy]) / d
        v = np.array([-u[1], u[0]])
        offsets = [(-1000.0, 0.0), (-800.0, 200.0), (-800.0, -200.0),
                   (-800.0, 2000.0), (-800.0, -2000.0),
                   (600.0, 200.0), (600.0, -200.0), (1000.0, 0.0)]
        for k, (a, b) in enumerate(offsets):
            mx, my = np.array([ox, oy]) + a * u + b * v
            x, y = cell(mx, my)
            lumlp[y, x] = 1000.0 + 100.0 * k
            altlp[y, x] = 10.0
    else:
        # A few lamps of one type near the centre (never on the observer pixel).
        lamps = [(c + 3, c + 2, 1.0), (c - 2, c + 4, 0.7), (c + 5, c - 3, 1.3),
                 (c - 4, c - 4, 0.9), (c + 1, c - 6, 1.1)]
        for x, y, flux in lamps:
            lumlp[y, x] = flux * 1000.0  # W/nm-ish; any positive value works
            altlp[y, x] = 10.0
        # One lamp on a tall mast to the north-east of the observer, placed so
        # that it lies inside the direct field of view (elevation 30 deg,
        # azimuth 45 deg). This makes the direct radiance/irradiance non-zero.
        dist = 4 * args.dx * np.sqrt(2.0)
        lumlp[c - 4, c + 4] = 500.0  # row c-4 is 4 px north (Fortran j = y_obs + 4)
        altlp[c - 4, c + 4] = 10.0 + dist * np.tan(np.deg2rad(30.0))

    # Extra lamps on a regular grid over the domain (--lamps N). The
    # observer pixel and the mast pixel are skipped. Deterministic.
    if args.lamps > 0:
        step = max(2, int(round(n / np.sqrt(args.lamps))))
        k = 0
        for y in range(step // 2, n, step):
            for x in range(step // 2, n, step):
                if (x, y) == (c, c) or (x, y) == (c + 4, c - 4):
                    continue
                if lumlp[y, x] == 0.0:
                    k += 1
                    lumlp[y, x] = 500.0 + 37.0 * (k % 11)
                    altlp[y, x] = 8.0 + (k % 5)

    if args.pad:
        if args.pad < n:
            sys.exit("--pad must be >= --size")
        before = (args.pad - n) // 2
        after = args.pad - n - before
        pad = lambda a: np.pad(a, ((before, after), (before, after)), "constant")
        topo, altlp, obsth, obstd, obstf, lumlp, origin = map(
            pad, (topo, altlp, obsth, obstd, obstf, lumlp, origin)
        )
        c += before
        n = args.pad

    b = args.basenm
    # Fortran index i (x) counts columns from 1, j (y) counts rows from the
    # bottom; twodin fills bindata(i, j) with j = nby..1 in file order.
    x_obs = c + 1
    y_obs = n - c
    save_bin(os.path.join(args.out, b + "_topogra.bin"), topo)
    save_bin(os.path.join(args.out, b + "_altlp.bin"), altlp)
    save_bin(os.path.join(args.out, b + "_obsth.bin"), obsth)
    save_bin(os.path.join(args.out, b + "_obstd.bin"), obstd)
    save_bin(os.path.join(args.out, b + "_obstf.bin"), obstf)
    save_bin(os.path.join(args.out, b + "_lumlp_001.bin"), lumlp)
    save_bin(os.path.join(args.out, "origin.bin"), origin)
    write_photometry(os.path.join(args.out, b + "_fctem_001.dat"))
    write_phase_file(os.path.join(args.out, "aerosol.txt"), ssa=0.9, g=0.7)
    write_phase_file(os.path.join(args.out, "layer.txt"), ssa=0.9, g=0.7)
    shutil.copy(
        os.path.join(REPO, "Molecular_optics", "MolecularAbs.txt"),
        os.path.join(args.out, "MolecularAbs.txt"),
    )
    write_illumina_in(os.path.join(args.out, "illumina.in"), b, args.dx, x_obs, y_obs, args)
    print("Case written to", args.out, "(%dx%d, dx=%g m, basenm=%s)" % (n, n, args.dx, b))
    if args.hill is not None:
        print("Hill: %g m, sigma %g m, centre (%g, %g) m; terrain at the observer %.1f m, max %.1f m"
              % (args.hill[0], args.hill[1], args.hill[2], args.hill[3], topo[c, c], topo.max()))
    print("Lamps: %d cells, observer cell (x=%d, y=%d)" % (np.count_nonzero(lumlp), x_obs, y_obs))


if __name__ == "__main__":
    main()
