"""Functions for analyzing an AIMD trajectory."""

import os
import sys
import glob
import numpy as np
import time
from sklearn.feature_selection import mutual_info_regression
from joblib import parallel_backend
import process


def charge_matrices() -> None:
    """
    Generates mutual information and corss-correlation matrices.

    """

    start_time = time.time()  # Used to report the executation speed
    # Search for the reference PDB
    pdbfile = process.get_pdb()

    # Set variables
    oldresi = "0"
    reslist = []
    atwbblist = []
    atwsclist = []

    # Search for the reference PDB
    with open(pdbfile, "r") as pdbfile:
        pdbfile_lines = pdbfile.readlines()

        for line in pdbfile_lines:
            if "ATOM" not in line:
                continue

            line_pieces = line.split()

            if line_pieces[4] != oldresi:
                reslist.append(line_pieces[3])
                oldresi = line_pieces[4]
                atwbblist.append([line_pieces[1]])
                atwsclist.append([])
                atom = line_pieces[2]

                if atom not in ["N", "O", "C", "H"]:
                    atwsclist.append([line_pieces[1]])
            else:
                atwbblist[int(oldresi) - 1].append(line_pieces[1])
                atom = line_pieces[2]

                if atom not in ["N", "O", "C", "H"]:
                    atwsclist[int(oldresi) - 1].append(line_pieces[1])

    with open("all_charges.xls", "r") as charge_file:
        # charge_file = open("all_charges.xls", "r").readlines()
        # Something must be breaking in this block
        charge_file_lines = charge_file.readlines()
        bbchargearray = []
        scchargearray = []
        charge_file_length = len(charge_file_lines)
        atwbblist_len = len(atwbblist)
        atwsclist_len = len(atwsclist)
        print(charge_file_length * 0.0005, "ps collected so far")

        for line2 in range(1, charge_file_length):
            charge_file_line_pieces = charge_file_lines[line2].split()
            bbchargearray.append([])
            scchargearray.append([])

            for resat in range(0, atwbblist_len):
                rescharge = 0.00

                for atwbbs in range(0, len(atwbblist[resat])):
                    rescharge += float(
                        charge_file_line_pieces[int(atwbblist[resat][atwbbs]) - 1]
                    )
                bbchargearray[line2 - 1].append(rescharge)

            for resat in range(0, atwsclist_len):
                ressccharge = 0.00

                for atwscs in range(0, len(atwsclist[resat])):
                    ressccharge += float(
                        charge_file_line_pieces[int(atwsclist[resat][atwscs]) - 1]
                    )
                scchargearray[line2 - 1].append(ressccharge)

    print("Finished looping through charges")

    np.array(bbchargearray)
    np.array(scchargearray)
    bbchargemutinf = np.array(bbchargearray).transpose()
    MI_mat = []
    bbchargearray_length = len(bbchargearray[0])

    print(f"Start looping through {bbchargearray_length}")

    for i in range(0, bbchargearray_length):
        print("rowMI", i)
        rowMI = mutual_info_regression(bbchargemutinf.transpose(), bbchargemutinf[i])
        MI_mat.append(rowMI)

    with open("mimatbb.csv", "w") as mimat:
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

    with open("chargematbb.csv", "w") as chargemat:
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

    with open("chargematsc.csv", "w") as chargemat:
        for resx in range(0, len(reslist)):
            for resy in range(0, len(reslist)):
                if resy == len(reslist) - 1:
                    extra = ""
                else:
                    extra = ","
                chargemat.write(str(corrcoefsc[resx][resy]) + extra)

            chargemat.write("\n")

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t-------------------------CHARGE MATRICES END--------------------------
        \tRESULT: Generated matrices.
        \tOUTPUT: Generated files.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )

    return


def get_joint_qres(res_x, res_y):
    """
    Calculates the joint net partial charge sum on each residue, q(RES).

    Parameters
    ----------
    res_x : str
        The amino acid to be represented on the x-axis
    res_y : str
        The amino acid to be represented on the y-axis

    Returns
    -------

    Notes
    -----
    Big picture flow of function.
        1. Read charges.xls in as a pd dataframe
        2. Get atom indices for requested atoms from get_res_indices()
        3. Extract and sum residue columns
        4. Save as a csv
        5. Return as a pandas dataframe
        6. Use joint_plot() to plot the results

    """

    start_time = time.time()  # Used to report the executation speed
    # Get the indices of the atoms to parse the charge.xls file
    residues = [res_x, res_y]
    for res in residues:
        atom_indices = process.get_res_atom_indices(res)
        with open("")




    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t-------------------------GET JOINT qRES END---------------------------
        \tRESULT: 
        \tOUTPUT: 
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )

if __name__ == "__main__":
    # Run when executed as a script
    get_joint_qres("Gln18", "Asp15")
