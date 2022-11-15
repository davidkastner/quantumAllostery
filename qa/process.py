"""Process the raw TeraChem output for future analysis."""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import time
from sklearn.feature_selection import mutual_info_regression


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
    start_time = time.time() # Used to report the executation speed
    frame_count = -1 # Report to user with -1 to account for header
    atoms = 489 # The number of atoms + 2 for header lines

    # Collect all qmscript.out files
    out_files = glob.glob("./**/qmscript.out", recursive=True)
    out_files_count = len(out_files) # Will report this to user
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
    
    total_time = round(time.time() - start_time, 3) # Seconds to run the function
    print(f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Combined {frame_count} frames from {out_files_count} runs.
        \tOUTPUT: Generated {all_charges} and {all_coors}.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """)


def charge_grab():
    """
    Placeholder function to show example docstring (NumPy format).

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from.

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution.
    """

    pdb = sys.argv[1]
    pdbfile = open(pdb + ".pdb", "r").readlines()
    chargefile = open("allcharge.xls", "r").readlines()
    oldresn = "NON"
    oldresi = "0"
    oldressc = "0"
    reslist = []
    atwbblist = []
    atwsclist = []
    atscolist = []

    for lines in range(0, len(pdbfile)):
        if "ATOM" in pdbfile[lines]:
            if pdbfile[lines].split()[5] != oldresi:
                reslist.append(pdbfile[lines].split()[3])
                oldresi = pdbfile[lines].split()[5]
                oldressc = pdbfile[lines].split()[5]
                atwbblist.append([pdbfile[lines].split()[1]])
                atwsclist.append([])
                if pdbfile[lines].split()[2] != "N":
                    if pdbfile[lines].split()[2] != "O":
                        if pdbfile[lines].split()[2] != "C":
                            if pdbfile[lines].split()[2] != "H":
                                atwsclist.append([pdbfile[lines].split()[1]])
            else:
                atwbblist[int(oldresi) - 1].append(pdbfile[lines].split()[1])
                if pdbfile[lines].split()[2] != "N":
                    if pdbfile[lines].split()[2] != "O":
                        if pdbfile[lines].split()[2] != "C":
                            if pdbfile[lines].split()[2] != "H":
                                atwsclist[int(oldresi) - 1].append(
                                    pdbfile[lines].split()[1]
                                )

    # Sidechains
    bbchargearray = []
    scchargearray = []
    print(len(chargefile) * 0.0005, "ps collected so far")
    for line2 in range(1, len(chargefile)):
        bbchargearray.append([])
        scchargearray.append([])
        for resat in range(0, len(atwbblist)):
            rescharge = 0.00
            for atwbbs in range(0, len(atwbblist[resat])):
                rescharge += float(
                    chargefile[line2].split()[int(atwbblist[resat][atwbbs]) - 1]
                )
            bbchargearray[line2 - 1].append(rescharge)
        for resat in range(0, len(atwsclist)):
            ressccharge = 0.00
            for atwscs in range(0, len(atwsclist[resat])):
                ressccharge += float(
                    chargefile[line2].split()[int(atwsclist[resat][atwscs]) - 1]
                )
            scchargearray[line2 - 1].append(ressccharge)

    np.array(bbchargearray)
    np.array(scchargearray)
    bbchargemutinf = np.array(bbchargearray).transpose()
    MI_mat = []

    for i in range(0, len(bbchargearray[0])):
        rowMI = mutual_info_regression(bbchargemutinf.transpose(), bbchargemutinf[i])
        MI_mat.append(rowMI)
    print(MI_mat)
    mimat = open("mimatbb.csv", "w")

    for resx in range(0, len(reslist)):
        for resy in range(0, len(reslist)):
            if resy == len(reslist) - 1:
                extra = ""
            else:
                extra = ","
            mimat.write(str(MI_mat[resx][resy]) + extra)
        mimat.write("\n")

    maxchgbb = np.amax(bbchargearray, axis=0)
    minchgbb = np.amin(bbchargearray, axis=0)
    avgchgbb = np.mean(bbchargearray, axis=0)
    stdchgbb = np.std(bbchargearray, axis=0)
    corrcoef = np.corrcoef(np.transpose(bbchargearray))
    maxchgsc = np.amax(scchargearray, axis=0)
    minchgsc = np.amin(scchargearray, axis=0)
    avgchgsc = np.mean(scchargearray, axis=0)
    stdchgsc = np.std(scchargearray, axis=0)
    corrcoefsc = np.corrcoef(np.transpose(scchargearray))
    print("BB stats below.")

    for resi in range(0, len(reslist)):
        print(
            f"{reslist[resi]} {maxchgbb[resi]} {minchgbb[resi]} {maxchgbb[resi] - minchgbb[resi]} {avgchgbb[resi]} {stdchgbb[resi]}"
        )

    chargemat = open("chargematbb.csv", "w")
    for resx in range(0, len(reslist)):
        for resy in range(0, len(reslist)):
            if resy == len(reslist) - 1:
                extra = ""
            else:
                extra = ","
            chargemat.write(str(corrcoef[resx][resy]) + extra)
        chargemat.write("\n")
    print("SC stats below.")

    for resi in range(0, len(reslist)):
        print(
            f"{reslist[resi]} {maxchgsc[resi]} {minchgsc[resi]} {maxchgsc[resi]-minchgsc[resi]} {avgchgsc[resi]} {stdchgsc[resi]}"
        )

    chargemat = open("chargematsc.csv", "w")
    for resx in range(0, len(reslist)):
        for resy in range(0, len(reslist)):
            if resy == len(reslist) - 1:
                extra = ""
            else:
                extra = ","
            chargemat.write(str(corrcoefsc[resx][resy]) + extra)
        chargemat.write("\n")

if __name__ == "__main__":
    # Run when executed as a script
    all_runs()