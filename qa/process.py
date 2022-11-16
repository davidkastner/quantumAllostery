"""Process the raw TeraChem output for future analysis."""

import os
import sys
import glob
import time


def all_runs():
    """
    Collects all charges or coordinates into single xls and xyz files.

    Likely the first executed function after generating the raw AIMD data.
    Trajectories were likely generated over multiple runs.
    This function combines all coordinate and charge data for each run.

    Notes
    -----
    Runs from inside the top-level directory of a TeraChem calculation.

    """

    # The files where all the charges and coors will be combined
    all_charges = f"all_charges.xls"
    all_coors = f"all_coors.xyz"
    start_time = time.time()  # Used to report the executation speed
    frame_count = -1  # Report to user with -1 to account for header
    atoms = 489  # The number of atoms + 2 for header lines

    # Collect all qmscript.out files
    out_files = glob.glob("./**/qmscript.out", recursive=True)
    out_files_count = len(out_files)  # Will report this to user
    run_info: list[list[int, str, str]] = []

    # Get the first MD step from out_files to identify where each job restarted
    # Restarting occurs from the last restart point not the last frame
    for out_file in out_files:
        out_content = open(out_file, "r").readlines()
        for line in out_content:
            if "MD STEP" in line:
                # A list of info about each run [start step, file, scr dir]
                md_step = int(line.split()[4])
                run_info.append([md_step, out_file])
                break

    # Get all src directories in safe, one for each restart attempt
    scrdir = glob.glob("./**/scr*/", recursive=True)
    scrdir.sort()  # Sort by modification date
    # Add the name of the scr directory location to the run_info list
    run_info.sort()
    for index, step in enumerate(run_info):
        step.append(scrdir[index])

    # Delete so we don't append to a previous version
    if os.path.exists(all_charges):
        os.remove(all_charges)
    if os.path.exists(all_coors):
        os.remove(all_coors)

    # Use the run_info to open the charge and coordinate files
    for index, run in enumerate(run_info):
        charge_file = open(f"{run[2]}charge.xls", "r").readlines()
        coors_file = open(f"{run[2]}coors.xyz", "r").readlines()
        # Combined charge and coors files
        all_charges_file = open(all_charges, "a")
        all_coors_file = open(all_coors, "a")

        # Parse runs differently based on if first, intermediate, or last
        if index == 0:
            for charge in range(0, run_info[index + 1][0] + 1):
                all_charges_file.write(charge_file[charge])
                frame_count += 1
            for coor in range(0, run_info[index + 1][0] * atoms):
                all_coors_file.write(coors_file[coor])
        elif (index < len(run_info) - 1) and index != 0:
            for charge in range(1, run_info[index + 1][0] - run[0] + 1):
                all_charges_file.write(charge_file[charge])
                frame_count += 1
            for coor in range(0, (run_info[index + 1][0] - run_info[index][0]) * atoms):
                all_coors_file.write(coors_file[coor])
        else:
            for charge in range(1, len(charge_file)):
                all_charges_file.write(charge_file[charge])
                frame_count += 1
            for coor in range(0, len(coors_file)):
                all_coors_file.write(coors_file[coor])

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Combined {frame_count} frames from {out_files_count} runs.
        \tOUTPUT: Generated {all_charges} and {all_coors}.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


if __name__ == "__main__":
    # Run when executed as a script
    all_runs()
