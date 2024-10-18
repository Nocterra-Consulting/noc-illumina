#!/usr/bin/env python3
#
# batch processing
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# March 2021

import os
import shutil
from collections import ChainMap, OrderedDict

from glob import glob
from itertools import product

import click
import numpy as np
import yaml


import illum
from illum import MultiScaleData as MSD
from illum.pytools import save_bin
import multiprocessing as mp
import tqdm as tqdm
import time as time

from functools import partial
from progressbar import progressbar
progress = partial(progressbar, redirect_stdout=True)


def input_line(val, comment, n_space=30):
    value_str = " ".join(str(v) for v in val)
    comment_str = " ; ".join(comment)
    return "%-*s ! %s" % (n_space, value_str, comment_str)


def MSDOpen(filename, cached={}):
    if filename in cached:
        return cached[filename]
    ds = MSD.Open(filename)
    cached[filename] = ds
    return ds


@click.command(name="batches")
@click.argument("input_path", type=click.Path(exists=True), default=".")
@click.argument("batch_name", required=False)
@click.option(
    "-c",
    "--compact",
    is_flag=True,
    help="If given, will chain similar executions. Reduces the overall number "
    "of runs at the cost of longuer individual executions.",
)
@click.option(
    "-N",
    "--batch_size",
    type=int,
    default=300,
    show_default=True,
    help="Number of runs per produced batch file.",
)
@click.option(
    "-s",
    "--scheduler",
    type=click.Choice(["parallel", "sequential", "slurm"]),
    default="sequential",
    help="Job scheduler",
)
def CLI_batches(input_path, compact, batch_size, scheduler, batch_name=None):
    """Makes the execution batches.

    INPUT_PATH is the path to the folder containing the inputs.

    BATCH_NAME is an optional name for the produced batch files.
    It overwrites the one defined in 'inputs_params.in' is given.
    """
    batches(input_path, compact, batch_size, scheduler, batch_name)


def batches(
    input_path=".",
    compact=False,
    batch_size=300,
    scheduler="sequential",
    batch_name=None,
):
    execute = dict(
        parallel="./execute &", sequential="./execute", slurm="sbatch ./execute"
    )
    execute_str = execute[scheduler] + "\n"

    os.chdir(input_path)

    with open("inputs_params.in") as f:
        params = yaml.safe_load(f)

    if batch_name is not None:
        params["batch_file_name"] = batch_name

    for fname in glob("%s*" % params["batch_file_name"]):
        os.remove(fname)

    exp_name = params["exp_name"]

    ds = MSD.Open(glob("*.hdf5")[0])

    # Pre process the obs extract
    print("Preprocessing...")
    shutil.rmtree("obs_data", True)
    lats, lons = ds.get_obs_pos()
    xs, ys = ds.get_obs_pos(proj=True)
    for lat, lon in zip(lats, lons):
        for i in range(len(ds)):
            os.makedirs("obs_data/%6f_%6f/%d" % (lat, lon, i))

    for i, fname in enumerate(progress(glob("*.hdf5")), 1):
        dataset = MSD.Open(fname)
        for clipped in dataset.split_observers():
            lat, lon = clipped.get_obs_pos()
            lat, lon = lat[0], lon[0]

            if "lumlp" in fname:
                clipped.set_buffer(0)
                clipped.set_overlap(0)
            for i, dat in enumerate(clipped):
                padded_dat = np.pad(dat, (512 - dat.shape[0]) // 2, "constant")
                save_bin(
                    "obs_data/%6f_%6f/%i/%s"
                    % (lat, lon, i, fname.rsplit(".", 1)[0] + ".bin"),
                    padded_dat,
                )
            if "srtm" in fname:
                for j in range(len(clipped)):
                    clipped[j][:] = 0
                clipped.save(f"obs_data/{lat:6f}_{lon:6f}/blank")

    # Add wavelength and multiscale
    spectral_bands = np.loadtxt("wav.lst", ndmin=2)
    params["wavelength"] = spectral_bands[:, 0].tolist()
    params["layer"] = list(range(len(ds)))
    params["observer_coordinates"] = list(zip(*ds.get_obs_pos()))

    wls = params["wavelength"]
    refls = np.loadtxt("refl.lst", ndmin=1).tolist()

    for pname in ["layer", "observer_coordinates"]:
        if len(params[pname]) == 1:
            params[pname] = params[pname][0]

    with open("lamps.lst") as f:
        lamps = f.read().split()

    if os.path.isfile("brng.lst"):
        brng = np.loadtxt("brng.lst", ndmin=1)
    else:
        brng = None

    # Clear and create execution folder
    dir_name = "exec" + os.sep
    shutil.rmtree(dir_name, True)
    os.makedirs(dir_name)

    count = 0
    multival = [k for k in params if isinstance(params[k], list)]
    multival = sorted(multival, key=len, reverse=True)  # Semi-arbitrary sort
    param_space = [params[k] for k in multival]
    N = np.prod([len(p) for p in param_space])

    
    # Create all possible combinations of parameters
    param_vals = list(product(*param_space))
    #args_list = list(zip(prod_params, func_array))
    #print(args_list)
    param_tasks = [OrderedDict(zip(multival, vals)) for vals in param_vals]
    const_args = [multival, brng, compact, dir_name, wls, 
                   refls, spectral_bands, lamps,  exp_name]
    const_str = ['multival', 'brng', 'compact', 'dir_name', 'wls', 
                 'refls', 'spectral_bands', 'lamps', 'exp_name']
    const_dict = dict(zip(const_str, const_args))
    # write to file yaml
    with open('const_params.yml', 'w') as f:
        yaml.dump(const_dict, f)
    # args_list = [(local_params, multival, params, brng, compact, dir_name, wls, 
    #                refls, spectral_bands, lamps,exp_name, ds) for local_params in param_tasks]
    time_one = time.time()
    # Setup folders
    folder_setup(param_vals, params, multival, brng, compact, dir_name, exp_name)
    time_two = time.time()
    print('Time to setup folders:', time_two-time_one)

    # write args_list to file
    print('Writing to file')
    with open('param_vals.txt', 'w') as f:
        for item in param_vals:
            f.write(f"{item}\n")
    with open('const_params.txt', 'w') as f:
            f.write(f"{const_args}")
    print('Done writing to file')
    
    return
    num_cpus = 20
    # Run multithreaded execution with a progress bar
    results = []
    with mp.Pool(processes=num_cpus) as pool:
        param_tasks = [OrderedDict(zip(multival, vals)) for vals in param_vals]
        args_list = [(local_params, multival, params, brng, compact, dir_name, wls, 
                   refls, spectral_bands, lamps, exe_location, exe_input, 
                   exe_output, exp_name, ds) for local_params in param_tasks]
        results = np.array(list(progress(
            pool.imap(execute_wrapper, args_list),
            total=N,
            desc="Processing"
        )))
        pool.close()
        pool.join()
    # for param_vals in progress(product(*param_space), max_value=N):
    #     param_generate(param_vals)
    ### NOCTERRA CHANGES   
    ### Write exes to text
    exe_input, exe_output, exe_location = results[0], results[1], results[2]
    with open('execute_info.txt', 'w') as f:
        for loc, input, output in zip(exe_location, exe_input, exe_output):
            input = input+'.in'
            output = output+'.out'
            f.write(f"{loc},{input},{output}\n")

    print("batches_testing")
    
    
    #print("Final count:", count)

    print("Done.")

def execute_wrapper(args: tuple) -> np.ndarray:
    """Wrapper to allow multiprocessing arguments to feed into execution."""
    return param_generate(*args)

def folder_setup(param_values, params, multival, brng, compact, dir_name, exp_name):
    print('folder setup')
    #NOCTERRA CHANGE
    ## Create array of locations and save names
    exe_location = []
    exe_input = []
    exe_output = []
    for param_vals in param_values:
        local_params = OrderedDict(zip(multival, param_vals))
        P = ChainMap(local_params, params)
        if (
            "azimuth_angle" in multival
            and P["elevation_angle"] == 90
            and params["azimuth_angle"].index(P["azimuth_angle"]) != 0
        ):
            continue

        if os.path.isfile("brng.lst"):
            obs_index = (
                0
                if "observer_coordinates" not in multival
                else params["observer_coordinates"].index(P["observer_coordinates"])
            )
            bearing = brng[obs_index]
        else:
            bearing = 0

        coords = "%6f_%6f" % P["observer_coordinates"]
        if "observer_coordinates" in multival:
            P["observer_coordinates"] = coords

        if compact:
            fold_name = (
                dir_name
                + os.sep.join(
                    f"{k}_{v}"
                    for k, v in local_params.items()
                    if k in ["observer_coordinates", "wavelength", "layer"]
                )
                + os.sep
            )
        else:
            fold_name = (
                dir_name
                + os.sep.join(f"{k}_{v}" for k, v in local_params.items())
                + os.sep
            )

        unique_ID = "-".join("%s_%s" % item for item in local_params.items())
        wavelength = "%g" % P["wavelength"]
        layer = P["layer"]
        reflectance = refls[wls.index(P["wavelength"])]
        bandwidth = spectral_bands[wls.index(P["wavelength"]), 1]
        create_symlinks(fold_name, params, exp_name, coords, layer, lamps, wavelength)

        ##NOCTERRA CHANGE
        exe_input.append(unique_ID)
        exe_output.append(f"{exp_name}_{unique_ID}")
        exe_location.append("%s" % os.path.abspath(fold_name))
        #print(unique_ID)
        ### NOCTERRA CHANGES
    print('finished loop')
    with open('execute_info.txt', 'w') as f:
        for loc, input, output in zip(exe_location, exe_input, exe_output):
            input = input+'.in'
            output = output+'.out'
            f.write(f"{loc},{input},{output}\n")
    return exe_location, exe_input, exe_output

def param_generate(local_params, multival, params, brng, compact, dir_name, wls, 
                   refls, spectral_bands, lamps, exe_location, exe_input, exe_output, 
                   exp_name, ds):
        P = ChainMap(local_params, params)
        wavelength = "%g" % P["wavelength"]
        layer = P["layer"]
        reflectance = refls[wls.index(P["wavelength"])]
        bandwidth = spectral_bands[wls.index(P["wavelength"]), 1]

        # Create symlinks
        create_symlinks(fold_name, params, exp_name, coords, layer, lamps, wavelength)
        # Create illumina.in
        input_data = create_illumina_in(exp_name, layer, ds, P, wavelength, reflectance, bandwidth, lamps, bearing )
        with open(fold_name + unique_ID + ".in", "w") as f:
            lines = (input_line(*zip(*line_data)) for line_data in input_data)
            f.write("\n".join(lines))


        # Write execute script
        ## NOCTERRA CHANGE Removing as we dont need
        # if not os.path.isfile(fold_name + "execute"):
        #     with open(fold_name + "execute", "w") as f:
        #         f.write("#!/bin/sh\n")
        #         f.write("#SBATCH --job-name=Illumina\n")
        #         f.write(
        #             "#SBATCH --time=%d:00:00\n" % params["estimated_computing_time"]
        #         )
        #         f.write("#SBATCH --mem=2G\n")
        #         f.write("cd %s\n" % os.path.abspath(fold_name))
        #         f.write("umask 0011\n")
        #     os.chmod(fold_name + "execute", 0o777)


            ## NOCTERRA CHANGE - DONT NEED TO CREATE BATCH LIST
            # Append execution to batch list
            # with open(f"{params['batch_file_name']}_{(count//batch_size)+1}", "a") as f:
            #     f.write("cd %s\n" % os.path.abspath(fold_name))
            #     f.write(execute_str)
            #     f.write("sleep 0.05\n")
            
            # count += 1

        # Add current parameters execution to execution script
        with open(fold_name + "execute", "a") as f:
            f.write("cp %s.in illumina.in\n" % unique_ID)
            f.write("./illumina\n")
            f.write(f"mv {exp_name}.out {exp_name}_{unique_ID}.out\n")
            f.write(f"mv {exp_name}_pcl.bin {exp_name}_pcl_{unique_ID}.bin\n")
    # return [exe_input, exe_output, exe_location]

def create_illumina_in(exp_name, layer, ds, P, wavelength, reflectance, bandwidth, lamps, bearing):
    # Create illumina.in
    input_data = (
        (("", "Input file for ILLUMINA"),),
        ((exp_name, "Root file name"),),
        (
            (ds.pixel_size(layer), "Cell size along X [m]"),
            (ds.pixel_size(layer), "Cell size along Y [m]"),
        ),
        (("aerosol.txt", "Aerosol optical cross section file"),),
        (
            ("layer.txt", "Layer optical cross section file"),
            (P["layer_aod"], "Layer aerosol optical depth at 500nm"),
            (P["layer_alpha"], "Layer angstom coefficient"),
            (P["layer_height"], "Layer scale height [m]"),
        ),
        ((P["double_scattering"] * 1, "Double scattering activated"),),
        ((P["single_scattering"] * 1, "Single scattering activated"),),
        ((wavelength, "Wavelength [nm]"), (bandwidth, "Bandwidth [nm]")),
        ((reflectance, "Reflectance"),),
        ((P["air_pressure"], "Ground level pressure [kPa]"),),
        (
            (P["aerosol_optical_depth"], "Aerosol optical depth at 500nm"),
            (P["angstrom_coefficient"], "Angstrom exponent"),
            (P["aerosol_height"], "Aerosol scale height [m]"),
        ),
        ((len(lamps), "Number of source types"),),
        ((P["stop_limit"], "Contribution threshold"),),
        (("", ""),),
        (
            (256, "Observer X position"),
            (256, "Observer Y position"),
            (P["observer_elevation"], "Observer elevation above ground [m]"),
        ),
        ((P["observer_obstacles"] * 1, "Obstacles around observer"),),
        (
            (P["elevation_angle"], "Elevation viewing angle"),
            ((P["azimuth_angle"] + bearing) % 360, "Azimuthal viewing angle"),
        ),
        ((P["direct_fov"], "Direct field of view"),),
        (("", ""),),
        (("", ""),),
        (("", ""),),
        (
            (
                P["reflection_radius"],
                "Radius around light sources where reflextions are computed",
            ),
        ),
        (
            (
                P["cloud_model"],
                "Cloud model: "
                "0=clear, "
                "1=Thin Cirrus/Cirrostratus, "
                "2=Thick Cirrus/Cirrostratus, "
                "3=Altostratus/Altocumulus, "
                "4=Cumulus/Cumulonimbus, "
                "5=Stratocumulus",
            ),
            (P["cloud_base"], "Cloud base altitude [m]"),
            (P["cloud_fraction"], "Cloud fraction"),
        ),
        (("", ""),),
    )
    return input_data

def create_symlinks(fold_name, params, exp_name, coords, layer, lamps, wavelength):
    if not os.path.isdir(fold_name):
            os.makedirs(fold_name)
            # Linking files
            mie_file = "{}_{}.txt".format(params["aerosol_profile"], wavelength)
            os.symlink(os.path.relpath(mie_file, fold_name), fold_name + "aerosol.txt")
            layer_file = "{}_{}.txt".format(params["layer_type"], wavelength)
            os.symlink(os.path.relpath(layer_file, fold_name), fold_name + "layer.txt")

            os.symlink(
                os.path.relpath("MolecularAbs.txt", fold_name),
                fold_name + "MolecularAbs.txt",
            )

            for i, lamp in enumerate(lamps, 1):
                os.symlink(
                    os.path.relpath(
                        f"fctem_wl_{wavelength}_lamp_{lamp}.dat",
                        fold_name,
                    ),
                    fold_name + exp_name + "_fctem_%03d.dat" % i,
                )

            illumpath = os.path.dirname(illum.__path__[0])
            os.symlink(
                os.path.abspath(illumpath + "/bin/illumina"), fold_name + "illumina"
            )

            # Copying layer data
            obs_fold = os.path.join("obs_data", coords, str(layer))

            os.symlink(
                os.path.relpath(os.path.join(obs_fold, "srtm.bin"), fold_name),
                fold_name + exp_name + "_topogra.bin",
            )

            os.symlink(
                os.path.relpath(os.path.join(obs_fold, "origin.bin"), fold_name),
                fold_name + "origin.bin",
            )

            for name in ["obstd", "obsth", "obstf", "altlp"]:
                os.symlink(
                    os.path.relpath(
                        os.path.join(obs_fold, f"{exp_name}_{name}.bin"),
                        fold_name,
                    ),
                    fold_name + f"{exp_name}_{name}.bin",
                )

            for i, lamp in enumerate(lamps, 1):
                os.symlink(
                    os.path.relpath(
                        os.path.join(
                            obs_fold,
                            f"{exp_name}_{wavelength}_lumlp_{lamp}.bin",
                        ),
                        fold_name,
                    ),
                    fold_name + "%s_lumlp_%03d.bin" % (exp_name, i),
                )
        