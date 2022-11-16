"""Analyze a AIMD trajectory."""

import os
import sys
import glob
import numpy as np
import time
from sklearn.feature_selection import mutual_info_regression


def charge_matrices():
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

    # Search for the reference PDB
    start_time = time.time() # Used to report the executation speed
    pdbfiles = glob.glob("./**/*.pdb", recursive=True)
    if len(pdbfiles) == 1:
        pdbfile = open(pdbfiles[0], "r").readlines()
    elif pdbfiles > 1:
        sys.exit('More the one PDB file was found.')
    else:
        sys.exit('No PDB files was found.')
    
    # Set variables
    charge_file = open("all_charges.xls", "r").readlines()
    oldresn = "NON"
    oldresi = "0"
    oldressc = "0"
    reslist = []
    atwbblist = []
    atwsclist = []
    atscolist = []

    for index,line in enumerate(pdbfile):
        if "ATOM" in line:
            if line.split()[5] != oldresi:
                reslist.append(line.split()[3])
                oldresi = line.split()[5]
                oldressc = line.split()[5]
                atwbblist.append(line.split()[1])
                atwsclist.append([])
                if line.split()[2] != "N":
                    if line.split()[2] != "O":
                        if line.split()[2] != "C":
                            if line.split()[2] != "H":
                                atwsclist.append(line.split()[1])
            else:
                atwbblist[int(oldresi) - 1].append(line.split()[1])
                if line.split()[2] != "N":
                    if line.split()[2] != "O":
                        if line.split()[2] != "C":
                            if line.split()[2] != "H":
                                atwsclist[int(oldresi) - 1].append(
                                    line.split()[1]
                                )

    # Sidechains
    bbchargearray = []
    scchargearray = []
    print(len(charge_file) * 0.0005, "ps collected so far")
    for line2 in range(1, len(charge_file)):
        bbchargearray.append([])
        scchargearray.append([])
        for resat in range(0, len(atwbblist)):
            rescharge = 0.00
            for atwbbs in range(0, len(atwbblist[resat])):
                rescharge += float(
                    charge_file[line2].split()[int(atwbblist[resat][atwbbs]) - 1]
                )
            bbchargearray[line2 - 1].append(rescharge)
        for resat in range(0, len(atwsclist)):
            ressccharge = 0.00
            for atwscs in range(0, len(atwsclist[resat])):
                ressccharge += float(
                    charge_file[line2].split()[int(atwsclist[resat][atwscs]) - 1]
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

    total_time = round(time.time() - start_time, 3) # Seconds to run the function
    print(f"""
        \t-------------------------CHARGE MATRICES END--------------------------
        \tRESULT: Generated matrices.
        \tOUTPUT: Generated files.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """)

if __name__ == "__main__":
    # Run when executed as a script
   charge_grab()