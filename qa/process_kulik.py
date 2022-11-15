"""Process the raw TeraChem output for future analysis."""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.feature_selection import mutual_info_regression


def xyz_consolidator():
    """
    Collects all charges and coordinates into single csv and xyz files.

    """

    filelist = glob.glob("*out*")
    print(filelist)
    fileinfo = []
    for files in filelist:
        myfile = open(files, "r").readlines()
        counter = 0
        for line in range(0, len(myfile)):
            if "MD STEP" in myfile[line]:
                if counter == 0:
                    fileinfo.append([int(myfile[line].split()[4]), files])
                    counter = 1

    fileinfosort = sorted(fileinfo)
    scrdir = glob.glob("safe/*/")
    scrdir.sort(key=os.path.getmtime)

    for lines in range(0, len(fileinfosort)):
        fileinfosort[lines].append(scrdir[lines])
    nat = 306

    for lines in range(0, len(fileinfosort)):
        chargefile = open(fileinfosort[lines][2] + "charge.xls", "r").readlines()
        coorsfile = open(fileinfosort[lines][2] + "coors.xyz", "r").readlines()
        allcharge = open("safe/allcharge.xls", "a")
        allcoors = open("safe/allcoors.xyz", "a")

        if lines == 0:
            for line2 in range(0, fileinfosort[lines + 1][0] + 1):
                allcharge.write(chargefile[line2])
            for line3 in range(0, fileinfosort[lines + 1][0] * nat):
                allcoors.write(coorsfile[line3])

        if lines < len(fileinfosort) - 1:
            if lines != 0:
                for line2 in range(
                    1, fileinfosort[lines + 1][0] - fileinfosort[lines][0] + 1
                ):
                    allcharge.write(chargefile[line2])
                for line3 in range(
                    0, (fileinfosort[lines + 1][0] - fileinfosort[lines][0]) * nat
                ):
                    allcoors.write(coorsfile[line3])

        if lines == len(fileinfosort) - 1:
            for line2 in range(1, len(chargefile)):
                allcharge.write(chargefile[line2])
            for line3 in range(0, len(coorsfile)):
                allcoors.write(coorsfile[line3])


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
