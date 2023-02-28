"""Manages file manipulations."""

import os
import time
import shutil
import glob


def run_all_replicates(function) -> None:
    """
    Loops through all replicates and runs a specific funtion on each.

    Often a function needs to be run on all replicates.
    Given a common file structure, this will take care of all of them at once.

    Notes
    -----
    The directories in the first level should contain replicates,
    and the sub directories should contain restarts.
    Additional random directories will lead to errors.
    """

    start_time = time.time()  # Used to report the executation speed

    # Set path elements
    root = os.getcwd()
    dirs = sorted(os.listdir(root))
    replicates_count = 0

    # Ignore directories
    ignore = ["Analysis"]

    # Loop over the replicate directories
    for dir in dirs:
        if dir in ignore:
            continue
        else:
            new_dir = f"{root}/{dir}"
            if os.path.isdir(new_dir):
                os.chdir(new_dir)
                print(f"> Executing function in {new_dir}.")
                function()
                replicates_count += 1

    total_time = round(time.time() - start_time, 3)  # Seconds to run
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Combined restarts for {replicates_count} replicates.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


def find_stalled():
    """
    Identifies frozen TeraChem jobs.

    When running TeraChem on SuperCloud,
    the jobs will occassionally freeze on the COSMO step.
    These instances can be difficult to identify.
    This script will check all jobs,
    and find those on the specific problematic step.
    """

    # Get the directories of each replicate
    primary = os.getcwd()
    dirs = sorted(glob.glob("*/"))
    stalled_jobs = []  # List of jobs currently on the COSMO step
    ignore = ["Analyze"]

    for dir in dirs:
        if dir in ignore:
            continue
        else:
            os.chdir(primary)  # Move back to replicates directory
            os.chdir(dir)
            secondary = os.getcwd()
            # Get the directories of each TeraChem frame calculation
            sub_dirs = sorted(glob.glob("*/"))

            for sub_dir in sub_dirs:
                if sub_dir in ignore:
                    continue
                else:
                    os.chdir(secondary)  # Move back to frames directory
                    os.chdir(sub_dir)
                    out_name = glob.glob("*.out")

                    if len(out_name) > 0:  # No out file if the job hasn't run
                        # Open the .out file
                        with open(out_name[0], "r") as out_file:
                            last_line = out_file.readlines()[-1]
                            if "=COSMO= size of" in last_line:
                                stalled_jobs.append(f"{dir}{sub_dir}")

    # Report findings to user
    if len(stalled_jobs) > 0:
        for index, stalled_job in enumerate(stalled_jobs):
            print(f"> Job in {stalled_jobs[index]} is likely stalled.")
    else:
        print("> No stalled jobs.")


def check_folder(dir_name) -> None:
    """
    Checks if a folder exists and creates it if it doesn't.

    Parameters
    ----------
    dir_name: str
        The name of the folder you are checking for
    """

    # Check if a folder exists, if it doesn't create one
    isExist = os.path.exists(dir_name)
    if not isExist:
        print(f"> The folder {dir_name} does not exist.")
        print(f"> Creating an {dir_name} folder")
        os.makedirs(dir_name)


def check_file(file_name, location) -> None:
    """
    Checks if a file exists and copies it if it doesn't.

    Parameters
    ----------
    file_name: str
        The name of the file you are checking for
    """

    # Check if a file exists, if it doesn't create one
    isExist = os.path.isfile(file_name)
    if not isExist:
        print(f"> Can't find {file_name}.")
        print(f"> Copying it over.")
        shutil.copyfile(location, f"{os.getcwd()}/{file_name}")


def copy_script(script_name) -> None:
    """
    Copies a script from the Scripts folder in the package.

    To perform a number of different analyses,
    you will require one of the saved scripts.
    This function will copy a requested script to your current location.

    Parameters
    ----------
    script_name: str
        The name of the requested script.

    """

    # The location of the cli
    this_dir = os.path.dirname(__file__)
    script_loc = f"{this_dir}/scripts/{script_name}"
    destination = f"{os.getcwd()}/{script_name}"
    shutil.copyfile(script_loc, destination)
