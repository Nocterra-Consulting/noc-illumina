#!/usr/bin/env python3

"""Tests for the photometry readers.

The unit tests build small synthetic files, so they run without the client
photometry library. The corpus tests read the real client library and skip
when that library is absent.
"""

import os

import numpy as np
import pytest

from illum import AngularPowerDistribution as APD

# ---------------------------------------------------------------- fixtures

# An isotropic source of 1000 lm emits 1000 / (4 pi) candela in every
# direction, so the flux integral of every fixture below is 1000 lm.
ISOTROPIC = 1000.0 / (4.0 * np.pi)

IES_VERTICAL = np.arange(0, 181, 10.0)


def write_ies(
    path,
    /,
    *,
    tilt="NONE",
    ballast=1.0,
    ballast_lamp=1.0,
    multiplier=1.0,
    keyword="[TEST] plain",
    encoding="utf-8",
    type=1,
    lumens=1000.0,
):
    """Write a synthetic LM-63 file that holds an isotropic distribution."""
    value = ISOTROPIC / (multiplier * ballast * ballast_lamp)
    lines = [
        "IESNA:LM-63-2002",
        keyword,
        f"TILT={tilt}",
    ]
    if tilt == "INCLUDE_1":
        lines[-1] = "TILT=INCLUDE"
        lines += ["1", "1", "0.0", "1.0"]
    elif tilt == "INCLUDE_13":
        lines[-1] = "TILT=INCLUDE"
        angles = " ".join(f"{a:.1f}" for a in np.linspace(0, 90, 13))
        factors = " ".join("1.0" for _ in range(13))
        lines += ["2", "13", angles, factors]
    nV = len(IES_VERTICAL)
    lines += [
        f"1 {lumens} {multiplier} {nV} 1 {type} 1 0.0 0.0 0.0",
        f"{ballast} {ballast_lamp} 10.0",
        " ".join(f"{a:.1f}" for a in IES_VERTICAL),
        "0.0",
        " ".join(f"{value:.6f}" for _ in range(nV)),
    ]
    with open(path, "wb") as f:
        f.write(("\n".join(lines) + "\n").encode(encoding))


def write_ldt(path, /, *, isym, Mc, Dc, n_planes, Ng=3, flux=1000.0, lorl=100.0):
    """Write a synthetic EULUMDAT file that holds an isotropic distribution."""
    header = ["Test company", "1", str(isym), str(Mc), str(Dc), str(Ng), "90"]
    header += ["text"] * 5
    header += ["0"] * 9
    header += ["100", str(lorl), "1", "0", "1"]
    lamp = ["1", "LED", str(flux), "4000", "80", "10"]
    ratios = ["0"] * 10

    C = [f"{i * Dc:g}" for i in range(Mc)]
    G = [f"{g:g}" for g in np.linspace(0, 180, Ng)]
    # The stored values are candela per 1000 lumens.
    body = [f"{ISOTROPIC * 1000.0 / flux:.6f}"] * (n_planes * Ng)

    lines = header + lamp + ratios + C + G + body
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def flux(apd):
    return float(np.sum(apd.to_type_c().vertical_profile(integrated=True)))


# -------------------------------------------------------------- unit tests


@pytest.mark.parametrize(
    "tilt", ["NONE", "INCLUDE_1", "INCLUDE_13", "tiltfile.tlt", "NONE "]
)
def test_tilt_forms(tmp_path, tilt):
    path = tmp_path / "lamp.ies"
    write_ies(path, tilt=tilt)
    apd = APD.from_ies(path)
    assert np.allclose(apd.vertical_angles, IES_VERTICAL)
    assert apd.lumens == pytest.approx(1000.0)
    assert flux(apd) == pytest.approx(1000.0, rel=1e-6)


def test_latin1_keyword(tmp_path):
    path = tmp_path / "latin1.ies"
    write_ies(path, keyword="[MORE] beam angle 120\xb0", encoding="latin-1")
    with pytest.raises(UnicodeDecodeError):
        open(path, encoding="utf-8").read()
    assert flux(APD.from_ies(path)) == pytest.approx(1000.0, rel=1e-6)


def test_ballast_factors_scale_the_candela(tmp_path):
    plain = tmp_path / "plain.ies"
    scaled = tmp_path / "scaled.ies"
    write_ies(plain)
    write_ies(scaled, ballast=0.9, ballast_lamp=1.1, multiplier=0.5)
    assert flux(APD.from_ies(plain)) == pytest.approx(
        flux(APD.from_ies(scaled)), rel=1e-5
    )


def test_float_lamp_count(tmp_path):
    path = tmp_path / "float.ies"
    write_ies(path)
    text = path.read_text().replace("\n1 1000.0 1.0", "\n1.0 1000.0 1.0")
    path.write_text(text)
    assert APD.from_ies(path).lumens == pytest.approx(1000.0)


def test_absolute_photometry_declares_no_flux(tmp_path):
    path = tmp_path / "absolute.ies"
    write_ies(path, lumens=-1)
    assert APD.from_ies(path).lumens == -1


def test_a_file_that_is_no_ies_raises_and_names_itself(tmp_path):
    path = tmp_path / "spectrum.ies"
    path.write_text("wavelength\trelativeIntensity\n330.00\t1.0\n331.00\t2.0\n")
    with pytest.raises(APD.PhotometryError, match="spectrum.ies"):
        APD.from_ies(path)


def test_type_a_raises_and_names_itself(tmp_path):
    path = tmp_path / "typea.ies"
    write_ies(path, type=3)
    with pytest.raises(APD.PhotometryError, match="typea.ies"):
        APD.from_ies(path)


def test_ldt_isym_4_with_more_c_angles_than_planes(tmp_path):
    path = tmp_path / "quadrant.ldt"
    write_ldt(path, isym=4, Mc=24, Dc=15, n_planes=7)
    apd = APD.from_ldt(path)
    assert apd.data.shape == (3, 7)
    assert np.allclose(apd.horizontal_angles, np.arange(0, 91, 15))
    assert apd.lumens == pytest.approx(1000.0)
    assert flux(apd) == pytest.approx(1000.0, rel=1e-5)


def test_ldt_that_declares_one_plane_but_stores_them_all(tmp_path):
    path = tmp_path / "allplanes.ldt"
    write_ldt(path, isym=1, Mc=24, Dc=15, n_planes=24)
    with pytest.warns(UserWarning, match="Isym 1"):
        apd = APD.from_ldt(path)
    # The reader closes the circle, so it adds the 360 degree plane.
    assert apd.data.shape == (3, 25)
    assert flux(apd) == pytest.approx(1000.0, rel=1e-5)


def test_ldt_light_output_ratio_scales_the_declared_flux(tmp_path):
    path = tmp_path / "lorl.ldt"
    write_ldt(path, isym=1, Mc=24, Dc=15, n_planes=1, lorl=80.0)
    assert APD.from_ldt(path).lumens == pytest.approx(800.0)


def test_a_file_that_is_no_ldt_raises_and_names_itself(tmp_path):
    path = tmp_path / "broken.ldt"
    path.write_text("\n".join(["text"] * 60) + "\n")
    with pytest.raises(APD.PhotometryError, match="broken.ldt"):
        APD.from_ldt(path)


def test_read_dispatches_on_the_extension(tmp_path):
    ies = tmp_path / "lamp.IES"
    write_ies(ies)
    assert flux(APD.read(ies)) == pytest.approx(1000.0, rel=1e-6)

    ldt = tmp_path / "lamp.LDT"
    write_ldt(ldt, isym=1, Mc=24, Dc=15, n_planes=1)
    assert flux(APD.read(ldt)) == pytest.approx(1000.0, rel=1e-5)

    with pytest.raises(APD.PhotometryError, match="lamp.xyz"):
        APD.read(tmp_path / "lamp.xyz")


# ------------------------------------------------------------ text round trip


def make_profile():
    """A downward beam that carries a small amount of uplight."""
    va = np.arange(181.0)
    data = 1000.0 * np.exp(-((va / 30.0) ** 2)) + 1e-6
    return APD.AngularPowerDistribution(
        vertical_angles=va,
        horizontal_angles=np.array([0.0]),
        data=data[:, None],
    )


def test_to_txt_writes_the_zenith_angle_in_column_two(tmp_path):
    path = tmp_path / "beam.lop"
    apd = make_profile()
    apd.to_txt(path)
    value, angle = np.loadtxt(path).T
    assert angle[0] == 0 and angle[-1] == 180
    # Angle 0 is the zenith, which is vertical angle 180 from the nadir.
    assert value[0] == pytest.approx(apd.data[180, 0])
    assert value[-1] == pytest.approx(apd.data[0, 0])


def test_txt_round_trip_keeps_the_sign_and_the_values(tmp_path):
    path = tmp_path / "beam.lop"
    apd = make_profile()
    apd.to_txt(path)
    back = APD.from_txt(path)

    assert np.all(np.diff(back.vertical_angles) > 0)
    assert np.allclose(back.vertical_angles, apd.vertical_angles)
    assert np.allclose(back.data, apd.data, rtol=1e-8)
    assert np.all(back.data >= 0)
    assert np.sum(back.vertical_profile(integrated=True)) > 0
    assert np.all(back.normalize().data >= 0)
    assert np.sum(back.normalize().vertical_profile(integrated=True)) == pytest.approx(
        1.0
    )


def test_from_txt_accepts_a_descending_angle_column(tmp_path):
    path = tmp_path / "descending.lop"
    apd = make_profile()
    apd.to_txt(path)
    rows = np.loadtxt(path)[::-1]
    np.savetxt(path, rows, fmt="%.9e")
    back = APD.from_txt(path)
    assert np.all(np.diff(back.vertical_angles) > 0)
    assert np.allclose(back.data, apd.data, rtol=1e-8)


def test_the_default_format_keeps_small_uplight(tmp_path):
    # The old notebooks wrote a peak-normalised profile with two decimals.
    apd = make_profile()
    apd.data = apd.data / apd.data.max()
    full = tmp_path / "full.lop"
    coarse = tmp_path / "coarse.lop"
    apd.to_txt(full)
    apd.to_txt(coarse, fmt="%.2f")

    def uplight(path):
        back = APD.from_txt(path)
        integral = back.vertical_profile(integrated=True)
        return np.sum(integral[back.vertical_angles > 90]) / np.sum(integral)

    assert uplight(coarse) == 0.0
    assert uplight(full) > 0.0


# ------------------------------------------------------------- corpus tests

GRIFFIN = (
    "/mnt/c/Nocterra OneDrive/OneDrive - Nocterra/01 Projects/N043 28South"
    "/N04301 Griffin Sports Complex Lighting Assessment/03 Project Info/ies"
)
ONEDRIVE = "/mnt/c/Nocterra OneDrive"


def corpus(root):
    """Collect every photometry file below a directory in one walk."""
    if not os.path.isdir(root):
        return []
    found = []
    for folder, _, names in os.walk(root):
        for name in names:
            if os.path.splitext(name)[1].lower() in (".ies", ".ldt"):
                found.append(os.path.join(folder, name))
    return found


# The walk over the client library is slow, so let a caller switch it off.
CORPUS = (
    []
    if os.environ.get("ILLUM_SKIP_CORPUS")
    else sorted(set(corpus(GRIFFIN) + corpus(ONEDRIVE)))
)

needs_corpus = pytest.mark.skipif(
    not CORPUS, reason="the client photometry library is absent"
)


@needs_corpus
@pytest.mark.parametrize("filename", CORPUS, ids=os.path.basename)
def test_every_corpus_file_parses_or_names_itself(filename):
    try:
        apd = APD.read(filename)
    except APD.PhotometryError as err:
        assert os.path.basename(filename) in str(err)
        return

    total = flux(apd)
    assert np.isfinite(total) and total > 0
    if apd.lumens > 0:
        assert total == pytest.approx(apd.lumens, rel=0.01)


@pytest.mark.skipif(
    not os.path.isfile(os.path.join(GRIFFIN, "vfl530.ies")),
    reason="the Griffin photometry set is absent",
)
def test_the_same_luminaire_reads_the_same_from_ies_and_ldt():
    ies = APD.read(os.path.join(GRIFFIN, "vfl530.ies"))
    ldt = APD.read(os.path.join(GRIFFIN, "vfl530.ldt"))

    angles = np.arange(181.0)
    a = np.interp(angles, ies.vertical_angles, ies.vertical_profile())
    b = np.interp(angles, ldt.vertical_angles, ldt.vertical_profile())

    peak = a > 0.01 * a.max()
    assert np.max(np.abs(b[peak] - a[peak]) / a[peak]) < 0.01
    assert np.sum(b) == pytest.approx(np.sum(a), rel=0.01)
