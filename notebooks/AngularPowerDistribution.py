#!/usr/bin/env python3

import os
import warnings
from dataclasses import dataclass

import matplotlib as mpl
import numpy as np
import scipy.interpolate


def mids(arr, /):
    return np.concatenate([[arr[0]], np.mean([arr[1:], arr[:-1]], 0), [arr[-1]]])


@dataclass
class AngularPowerDistribution:
    vertical_angles: np.ndarray
    horizontal_angles: np.ndarray
    data: np.ndarray
    type: int = 1
    lumens: float = -1
    name: str = ""

    def __repr__(self):
        return f"APD<{self.data.shape}:{self._type_letter()}>"

    def _type_letter(self):
        return "CBA"[self.type - 1]

    def cycle(self, *args, **kwargs):
        return cycle(self, *args, **kwargs)

    def interpolate(self, *args, **kwargs):
        return interpolate(self, *args, **kwargs)

    def normalize(self, *args, **kwargs):
        return normalize(self, *args, **kwargs)

    def plot1d(self, *args, **kwargs):
        return plot1d(self, *args, **kwargs)

    def plot2d(self, *args, **kwargs):
        return plot2d(self, *args, **kwargs)

    def plot3d(self, *args, **kwargs):
        return plot3d(self, *args, **kwargs)

    def to_type_c(self, *args, **kwargs):
        return to_type_c(self, *args, **kwargs)

    def to_ies(self, filename, *args, **kwargs):
        return to_ies(filename, self, *args, **kwargs)

    def to_txt(self, filename, *args, **kwargs):
        return to_txt(filename, self, *args, **kwargs)

    def vertical_profile(self, *args, **kwargs):
        return vertical_profile(self, *args, **kwargs)


class PhotometryError(ValueError):
    """A photometry file is unreadable or holds an unsupported format."""


def _read_text(filename, /):
    """Read a photometry file as text.

    Photometry files carry no encoding declaration. Try utf-8 first and fall
    back to latin-1, which accepts every byte value.
    """
    with open(filename, "rb") as f:
        raw = f.read()
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def from_ies(filename, /):
    """Read an IESNA LM-63 photometry file.

    The candela values scale by the candela multiplier, the ballast factor and
    the ballast-lamp photometric factor. Type C (1) and type B (2) files are
    accepted. Type A (3) raises a PhotometryError.

    A lamp flux of -1 marks absolute photometry. Such a file declares no flux,
    so `lumens` becomes -1.
    """
    lines = _read_text(filename).splitlines()
    for index, line in enumerate(lines):
        if line.upper().lstrip().startswith("TILT="):
            break
    else:
        raise PhotometryError(f"{filename}: the file holds no TILT= line, so it is no IES file.")

    tokens = " ".join(lines[index + 1 :]).replace(",", " ").split()

    if "INCLUDE" in line.upper():
        # The tilt block holds <geometry> <n> <n angles> <n multipliers>.
        try:
            n_tilt = int(float(tokens[1]))
        except (IndexError, ValueError):
            raise PhotometryError(f"{filename}: the TILT=INCLUDE block is unreadable.") from None
        if n_tilt < 0 or len(tokens) < 2 + 2 * n_tilt:
            raise PhotometryError(f"{filename}: the TILT=INCLUDE block is too short.")
        tokens = tokens[2 + 2 * n_tilt :]

    if len(tokens) < 13:
        raise PhotometryError(f"{filename}: the IES header is too short.")

    try:
        n_lamps = int(float(tokens[0]))
        lumens_per_lamp = float(tokens[1])
        multiplier = float(tokens[2])
        nV = int(float(tokens[3]))
        nH = int(float(tokens[4]))
        type = int(float(tokens[5]))
        ballast_factor = float(tokens[10])
        ballast_lamp_factor = float(tokens[11])
    except ValueError as err:
        raise PhotometryError(f"{filename}: the IES header is unreadable ({err}).") from None

    if nV < 1 or nH < 1:
        raise PhotometryError(f"{filename}: the IES header declares {nV} vertical and {nH} horizontal angles.")
    if type not in (1, 2, 3):
        raise PhotometryError(f"{filename}: the IES header declares photometric type {type}.")
    if type == 3:
        raise PhotometryError(f"{filename}: photometric type A is unsupported.")

    body = tokens[13:]
    if len(body) < nV + nH + nV * nH:
        raise PhotometryError(
            f"{filename}: the file holds {len(body)} values after the header, "
            f"but the header asks for {nV + nH + nV * nH}."
        )

    try:
        values = np.array(body[: nV + nH + nV * nH], dtype="float64")
    except ValueError as err:
        raise PhotometryError(f"{filename}: the candela block is unreadable ({err}).") from None

    factor = multiplier * ballast_factor * ballast_lamp_factor
    va = values[:nV]
    ha = values[nV : nV + nH]
    data = values[nV + nH :].reshape(nH, nV).T * factor

    # A lamp flux of -1 marks absolute photometry, so the file declares no flux.
    lumens = -1.0 if lumens_per_lamp < 0 else n_lamps * lumens_per_lamp

    return AngularPowerDistribution(
        lumens=lumens,
        type=type,
        vertical_angles=va,
        horizontal_angles=ha,
        data=data,
        name=str(filename),
    )


# The number of intensity planes that each EULUMDAT symmetry code implies.
def _ldt_planes_from_isym(isym, Mc):
    return {0: Mc, 1: 1, 2: Mc // 2 + 1, 3: Mc // 2 + 1, 4: Mc // 4 + 1}.get(isym)


def _ldt_split(filename, body, Mc, Ng, isym):
    """Find how many C angles and how many intensity planes the body holds.

    EULUMDAT files store `A` C angle values, `Ng` G angle values and then
    `P * Ng` intensity values. Real files disagree with their own Isym code,
    so derive A and P by search and validate the angle columns.
    """
    for A in (Mc, Mc // 4 + 1, Mc // 2 + 1, 1):
        if A < 1 or len(body) <= A + Ng:
            continue
        rest = len(body) - A - Ng
        if rest % Ng:
            continue
        P = rest // Ng
        if P < 1 or P > Mc or P > A:
            continue
        C = body[:A]
        G = body[A : A + Ng]
        if np.any(np.diff(C) <= 0) or C[0] < 0 or C[-1] >= 360:
            continue
        if np.any(np.diff(G) <= 0) or G[0] < -90 or G[-1] > 180:
            continue
        expected = _ldt_planes_from_isym(isym, Mc)
        if expected is not None and expected != P:
            warnings.warn(
                f"{filename}: Isym {isym} implies {expected} intensity planes, "
                f"but the file stores {P}.",
                stacklevel=3,
            )
        return A, P
    raise PhotometryError(
        f"{filename}: the EULUMDAT body holds {len(body)} values, which fits no "
        f"C angle count against Mc={Mc} and Ng={Ng}."
    )


def from_ldt(filename, /):
    """Read an EULUMDAT (.ldt) photometry file.

    EULUMDAT stores intensities in candela per 1000 lumens, so the values scale
    by the conversion factor of header line 24 and by `flux / 1000`. The flux
    comes from the first lamp set. A file that lists several lamp sets lists
    dimming steps, so the first set holds the full output. Some files mark a
    negative lamp count, which the standard reads as absolute photometry, but
    every such file in the test corpus still stores candela per 1000 lumens,
    so the flux scaling always applies.

    The reported `lumens` is the luminaire flux, that is the lamp flux times
    the light output ratio of header line 23. The intensity values already
    carry that ratio.

    The G angles are measured from the nadir. A file that stores G from -90 to
    90 is shifted by +90 to reach the 0 to 180 nadir convention.

    The azimuth origin is shifted so that the first stored C plane sits at 0.
    That rotation about the vertical axis leaves the vertical profile unchanged
    and it lets `cycle` mirror the stored planes correctly.
    """
    lines = _read_text(filename).splitlines()

    def number(index, cast=float):
        try:
            return cast(float(lines[index].strip().replace(",", ".")))
        except (IndexError, ValueError):
            raise PhotometryError(
                f"{filename}: line {index + 1} holds no number, so this is no EULUMDAT file."
            ) from None

    if len(lines) < 42:
        raise PhotometryError(f"{filename}: the file holds {len(lines)} lines, which is too few for EULUMDAT.")

    isym = number(2, int)
    Mc = number(3, int)
    Ng = number(5, int)
    conversion = number(23)
    n_sets = number(25, int)

    if Mc < 1 or Ng < 1 or n_sets < 1:
        raise PhotometryError(f"{filename}: the EULUMDAT header declares Mc={Mc}, Ng={Ng} and {n_sets} lamp sets.")

    # Each lamp set holds six lines: count, type, flux, colour temp, CRI, watts.
    # A file that lists several sets lists dimming steps, so the first set holds
    # the full output.
    flux = abs(number(26 + 2))
    lorl = number(22)

    start = 36 + 6 * n_sets
    try:
        body = np.array(" ".join(lines[start:]).replace(",", " ").split(), dtype="float64")
    except ValueError as err:
        raise PhotometryError(f"{filename}: the EULUMDAT body is unreadable ({err}).") from None

    A, P = _ldt_split(filename, body, Mc, Ng, isym)

    ha = body[:P].copy()
    va = body[A : A + Ng].copy()
    data = body[A + Ng : A + Ng + P * Ng].reshape(P, Ng).T * conversion

    if flux > 0:
        data = data * (flux / 1000.0)

    if va[0] < 0:
        if va[0] < -90 or va[-1] > 90:
            raise PhotometryError(f"{filename}: the G angles span {va[0]} to {va[-1]}, which fits no known convention.")
        va = va + 90.0

    ha = ha - ha[0]

    # Close the circle when the stored planes cover the full azimuth.
    if P > 1:
        step = ha[1] - ha[0]
        if np.isclose(ha[-1] + step, 360.0):
            ha = np.concatenate([ha, [360.0]])
            data = np.concatenate([data, data[:, :1]], axis=1)

    return AngularPowerDistribution(
        lumens=flux * lorl / 100.0,
        type=1,
        vertical_angles=va,
        horizontal_angles=ha,
        data=data,
        name=str(filename),
    )


def read(filename, /):
    """Read a photometry file and dispatch on its extension."""
    ext = os.path.splitext(str(filename))[1].lower()
    if ext == ".ies":
        return from_ies(filename)
    if ext == ".ldt":
        return from_ldt(filename)
    if ext in (".lop", ".txt", ".dat"):
        return from_txt(filename)
    raise PhotometryError(f"{filename}: the extension '{ext}' names no known photometry format.")


def to_type_c(apd, /, *, step=1, method="linear"):
    """Resample a type B distribution onto a type C grid.

    A type B sample at lateral angle B and vertical angle V points in the
    direction that satisfies `cos(gamma) = cos(B) * cos(V)`, where gamma is the
    polar angle from the nadir. Map every sample to a (gamma, azimuth) pair and
    resample onto the full type C grid.
    """
    if apd.type == 1:
        return apd
    if apd.type == 3:
        raise PhotometryError(f"{apd.name or 'this distribution'}: photometric type A is unsupported.")

    V = np.asarray(apd.vertical_angles, dtype="float64")
    B = np.asarray(apd.horizontal_angles, dtype="float64")
    data = np.asarray(apd.data, dtype="float64")

    # A file that stores only one side declares symmetry about that plane.
    if B[0] == 0 and len(B) > 1:
        B = np.concatenate([-B[:0:-1], B])
        data = np.concatenate([data[:, :0:-1], data], axis=1)
    if V[0] == 0 and len(V) > 1:
        V = np.concatenate([-V[:0:-1], V])
        data = np.concatenate([data[:0:-1], data], axis=0)

    Bg, Vg = np.meshgrid(np.deg2rad(B), np.deg2rad(V))
    gamma = np.rad2deg(np.arccos(np.clip(np.cos(Vg) * np.cos(Bg), -1, 1))).ravel()
    phi = np.rad2deg(np.arctan2(np.sin(Vg), -np.cos(Vg) * np.sin(Bg))).ravel() % 360
    values = data.ravel()

    N = round(90 / step)
    H = np.linspace(0, 360, 4 * N + 1)
    Vc = np.linspace(0, 180, 2 * N + 1)

    # The apex carries no azimuth, so repeat it across the whole grid.
    apex = gamma < 1e-9
    if np.any(apex):
        gamma = np.concatenate([gamma[~apex], np.zeros_like(H)])
        phi = np.concatenate([phi[~apex], H])
        values = np.concatenate([values[~apex], np.full_like(H, np.mean(values[apex]))])

    # Repeat the samples on both sides so the azimuth wraps.
    points_h = np.concatenate([phi - 360, phi, phi + 360])
    points_v = np.tile(gamma, 3)
    points_val = np.tile(values, 3)

    grid = scipy.interpolate.griddata(
        (points_h, points_v),
        points_val,
        tuple(np.meshgrid(H, Vc)),
        method=method,
        fill_value=0,
    )

    return AngularPowerDistribution(
        lumens=apd.lumens,
        type=1,
        vertical_angles=Vc,
        horizontal_angles=H,
        data=grid,
        name=apd.name,
    )


def to_ies(filename, apd, /):
    out = [
        "IES:LM-63-2019",
        "TILT=NONE",
        f"1 {apd.lumens} 1.0 {apd.data.shape[0]} {apd.data.shape[1]} "  # cont.
        f"{apd.type} 1 0 0 0",
        "1.0 1.0 0.0",
    ]
    with open(filename, "w") as f:
        f.write("\n".join(out) + "\n")
        apd.vertical_angles.tofile(f, sep=" ", format="%f")
        f.write("\n")
        apd.horizontal_angles.tofile(f, sep=" ", format="%f")
        f.write("\n")
        for row in apd.data.T:
            row.tofile(f, sep=" ", format="%f")
            f.write("\n")


def from_txt(filename, /):
    """Read a `.lop` profile.

    Column 1 holds the value and column 2 holds the angle from the zenith, so
    angle 0 points straight up and angle 180 points at the nadir. Rows may
    ascend or descend; the angles are sorted after the read.
    """
    data, ang = np.loadtxt(filename).T
    va = 180.0 - ang
    order = np.argsort(va, kind="stable")

    return AngularPowerDistribution(
        type=1,
        vertical_angles=va[order],
        horizontal_angles=np.array([0.0]),
        data=data[order][:, None],
        name=str(filename),
    )


def to_txt(filename, apd, /, *, fmt="%.9e", **kwargs):
    """Write a `.lop` profile.

    Column 1 holds the value and column 2 holds the angle from the zenith.
    The default format keeps full precision, because a coarse format rounds
    small uplight values to zero.
    """
    zenith = np.arange(181)
    data = np.interp(
        180 - zenith, apd.vertical_angles, apd.vertical_profile(), left=0, right=0
    )
    np.savetxt(filename, np.stack((data, zenith), axis=1), fmt=fmt, **kwargs)


def vertical_profile(apd, /, *, integrated=False):
    apd = to_type_c(apd)

    profile = (
        np.average(
            apd.data,
            axis=1,
            weights=np.diff(mids(apd.horizontal_angles)),
        )
        if len(apd.horizontal_angles) > 1
        else apd.data[:, 0].copy()
    )
    if integrated:
        profile *= 2 * np.pi * np.diff(-np.cos(np.deg2rad(mids(apd.vertical_angles))))
    return profile


def normalize(apd):
    data = apd.data / np.sum(apd.vertical_profile(integrated=True))

    return AngularPowerDistribution(
        type=apd.type,
        vertical_angles=apd.vertical_angles,
        horizontal_angles=apd.horizontal_angles,
        data=data,
        name=apd.name,
    )


def cycle(apd, /, *, step=1, kind="linear"):
    apd = to_type_c(apd)

    ha = apd.horizontal_angles
    data = apd.data
    if len(ha) == 1:
        ha = np.array([0, 90])
        data = np.repeat(data, 2, axis=1)
    if ha[-1] == 90:
        ha = np.concatenate([ha[:-1], 180 - ha[::-1]])
        data = np.concatenate([data[:, :-1], data[:, ::-1]], axis=1)
    if ha[-1] == 180:
        ha = np.concatenate([ha[:-1], 360 - ha[::-1]])
        data = np.concatenate([data[:, :-1], data[:, ::-1]], axis=1)

    apd_cycle = AngularPowerDistribution(
        lumens=apd.lumens,
        vertical_angles=apd.vertical_angles,
        horizontal_angles=ha,
        data=data,
        name=apd.name,
    )

    return apd_cycle


def interpolate(apd, /, *, step=1, method="linear"):
    N = round(90 / step)
    apd = apd.cycle()
    ha, va = np.meshgrid(apd.horizontal_angles, apd.vertical_angles)
    H = np.linspace(0, 360, 4 * N + 1)
    V = np.linspace(0, 180, 2 * N + 1)

    interp = scipy.interpolate.griddata(
        (ha.flatten(), va.flatten()),
        apd.data.flatten(),
        tuple(np.meshgrid(H, V)),
        method=method,
        fill_value=0,
    )

    return AngularPowerDistribution(
        lumens=apd.lumens,
        vertical_angles=V,
        horizontal_angles=H,
        data=interp,
        name=apd.name,
    )


def plot1d(apd, /, ax=None, **kwargs):
    if ax is None:
        ax = mpl.pyplot.gca()

    return ax.plot(apd.vertical_angles, apd.vertical_profile(), **kwargs)


def plot2d(
    apd,
    /,
    ax=None,
    *,
    wrap=False,
    interpolation="nearest",
    cmap=mpl.rc_params()["image.cmap"],
    vmin=0,
    vmax=None,
):
    ha = apd.horizontal_angles
    va = apd.vertical_angles
    data = apd.data

    if wrap and ha[-1] == 360:
        ha = np.roll(ha[:-1], len(ha) // 2)
        ha = (ha + 180) % 360 - 180
        ha = np.concatenate([ha, ha[0:1] + 360])
        data = np.roll(data[:, :-1], len(ha) // 2, axis=1)
        data = np.concatenate([data, data[:, 0:1]], axis=1)

    if ax is None:
        ax = mpl.pyplot.gca()

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    im = mpl.image.NonUniformImage(
        ax,
        cmap=cmap,
        norm=norm,
        interpolation=interpolation,
        extent=[ha[0], ha[-1], va[0], va[-1]],
    )

    im.set_data(ha, va, data)
    ax.add_image(im)
    ax.set_xlim(ha[0], ha[-1])
    ax.set_ylim(va[0], va[-1])

    return im


def plot3d(apd, /, *, wireframe=False, **kwargs):
    ha, va = np.meshgrid(
        np.deg2rad(90 - apd.horizontal_angles),
        np.deg2rad(180 - apd.vertical_angles),
        copy=False,
        sparse=True,
    )
    r = apd.data
    x = r * np.sin(va) * np.cos(ha)
    y = r * np.sin(va) * np.sin(ha)
    z = r * np.cos(va)

    m = max(np.max(np.abs(x)), np.max(np.abs(y)), np.max(np.abs(z)))

    fig = mpl.pyplot.figure()
    ax = fig.add_subplot(projection="3d")
    (ax.plot_wireframe if wireframe else ax.plot_surface)(
        x, y, z, rstride=1, cstride=1, **kwargs
    )
    ax.set_xlim(-m, m)
    ax.set_ylim(-m, m)
    ax.set_zlim(-m, m)
