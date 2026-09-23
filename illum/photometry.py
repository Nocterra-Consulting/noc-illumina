#!/usr/bin/env python3

import hashlib
import os
import shutil
import warnings

import click
import numpy as np

from illum import AngularPowerDistribution as APD


@click.command(name="photometry")
@click.argument("filenames", nargs=-1, type=click.Path(exists=True), required=True)
@click.option(
    "-o",
    "--outdir",
    default=".",
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option("--fmt", default="%.9e", help="Number format of the output values.")
def CLI_photometry(filenames, outdir, fmt):
    """Convert photometry files to the Illumina `.lop` format.

    FILENAMES lists the IES, EULUMDAT or `.lop` files to convert. The command
    writes one `.lop` file per input into the output folder and prints one
    summary line per file. The command exits with code 1 if any file fails.
    """
    photometry(filenames, outdir, fmt)


def photometry(filenames, outdir=".", fmt="%.9e"):
    os.makedirs(outdir, exist_ok=True)
    failed = 0

    for filename in filenames:
        name = os.path.basename(filename)
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                apd = APD.read(filename)
                messages = [str(w.message) for w in caught]

            grid = apd.to_type_c()
            profile = grid.vertical_profile(integrated=True)
            total = np.sum(profile)
            if not np.isfinite(total) or total <= 0:
                raise APD.PhotometryError(f"{filename}: the flux integral is {total}.")

            uplight = np.sum(profile[grid.vertical_angles > 90]) / total
            outname = os.path.join(outdir, os.path.splitext(name)[0] + ".lop")
            apd.to_txt(outname, fmt=fmt)
        except Exception as err:
            failed += 1
            click.echo(f"{name}: FAILED: {err}", err=True)
            continue

        ext = os.path.splitext(name)[1].lower().lstrip(".")
        click.echo(
            f"{name}: {ext} type {apd._type_letter()}, "
            f"{len(apd.horizontal_angles)} planes, uplight {uplight * 100:.3f}%"
        )
        for message in messages:
            click.echo(f"{name}: warning: {message}", err=True)

    if failed:
        click.echo(f"{failed} of {len(filenames)} files failed.", err=True)
        raise SystemExit(1)


NOTEBOOK = "photometry_inspection.ipynb"
READER = "AngularPowerDistribution.py"


@click.command(name="photometry-notebook")
@click.argument("outdir", type=click.Path(file_okay=False))
def CLI_photometry_notebook(outdir):
    """Write the standalone photometry inspection notebook.

    OUTDIR receives the notebook and a copy of the reader it needs. The pair
    then runs on a plain Python install, away from this checkout. Run this
    command again after every change to the reader, so the copy never drifts.
    """
    photometry_notebook(outdir)


def photometry_notebook(outdir):
    package = os.path.dirname(os.path.abspath(__file__))
    sources = {
        NOTEBOOK: os.path.join(os.path.dirname(package), "notebooks", NOTEBOOK),
        READER: os.path.join(package, READER),
    }
    for name, path in sources.items():
        if not os.path.isfile(path):
            raise click.ClickException(f"{path}: this copy of illum holds no {name}.")

    os.makedirs(outdir, exist_ok=True)
    for name, path in sources.items():
        shutil.copyfile(path, os.path.join(outdir, name))
        click.echo(f"wrote {os.path.join(outdir, name)}")

    with open(sources[READER], "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    click.echo(f"reader sha256 {digest}")
