"""Functions for analyzing an AIMD trajectory."""

import os
import sys
import glob
import shutil
import numpy as np
import time
from sklearn.feature_selection import mutual_info_regression
from joblib import parallel_backend
import pandas as pd
import qa.process
import qa.plot
import qa.manage


def charge_matrix() -> None:
    """
    Generates mutual information and cross-correlation matrices.

    """

    start_time = time.time()  # Used to report the executation speed
    # Search for the reference PDB
    pdbfile = qa.process.get_pdb()

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
        print(f"> {round(charge_file_length * 0.0005, 2)} ps collected so far")

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

    print("> Finished looping through charges")

    np.array(bbchargearray)
    np.array(scchargearray)
    bbchargemutinf = np.array(bbchargearray).transpose()
    MI_mat = []
    bbchargearray_length = len(bbchargearray[0])

    print(f"> Start looping through {bbchargearray_length} residues")

    for i in range(0, bbchargearray_length):
        print(f"   > rowMI {i}")
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

    # print("> BB stats below.")
    # for resi in range(0, len(reslist)):
    #     print(
    #         f"{reslist[resi]} {maxchgbb[resi]} {minchgbb[resi]} {maxchgbb[resi] - minchgbb[resi]} {avgchgbb[resi]} {stdchgbb[resi]}"
    #     )

    with open("chargematbb.csv", "w") as chargemat:
        for resx in range(0, len(reslist)):
            for resy in range(0, len(reslist)):
                if resy == len(reslist) - 1:
                    extra = ""
                else:
                    extra = ","
                chargemat.write(str(corrcoef[resx][resy]) + extra)

            chargemat.write("\n")

    # print("SC stats below.")
    # for resi in range(0, len(reslist)):
    #     print(
    #         f"{reslist[resi]} {maxchgsc[resi]} {minchgsc[resi]} {maxchgsc[resi]-minchgsc[resi]} {avgchgsc[resi]} {stdchgsc[resi]}"
    #     )

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
    # Add folder for the analysis
    root = os.getcwd()
    os.chdir("Analysis")
    if not os.path.isdir("3_coupling"):
        os.makedirs("3_coupling")
    os.chdir(root)

    # Create a new data frame to append the two residues of interest
    residues = [res_x, res_y]
    joint_df = pd.DataFrame(columns=residues)

    # Get the indices of the atoms to parse the charge.xls file
    for res in residues:
        # Search for a .xls file and read it in
        charge_file = qa.process.get_charge_file()
        charge_df = pd.read_csv(charge_file, sep="\t")

        # Get the atom indices of the residue
        atom_indices = qa.process.get_res_atom_indices(res)
        # Sum the charges of the atoms of the requested residues
        summed_charges = charge_df[charge_df.columns[atom_indices]].sum(axis=1).tolist()
        charge_df.to_csv("charge_df.csv")
        
        # Add the summed residue to the new dataframe
        joint_df[res] = summed_charges

    ext = "png"
    plot_name = f"{res_x}_{res_y}_dist.{ext}"
    qa.plot.get_charge_distributions(joint_df, plot_name, res_x, res_y, ext)

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t-------------------------GET JOINT qRES END-------------------------
        \tRESULT: Extracted and computed the joint charge distribution.
        \tOUTPUT: Created the charge distribution plot: {plot_name}.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )

    return joint_df

def cpptraj_covars(delete, recompute=False) -> None:
    """
    Calculate the covariance using CPPTraj.

    A measure of how different amino acids are geometrically correlated.
    For example, fluctuations in side chain position or conformational changes.

    Parameters
    ----------
    delete: List[List[int]]
        A list containing two lists. The first is the residues to delete.
        The second is the columns of the matrix to delete.
    recompute: bool
        Recompute the calculation even if results are already present.

    """
    job = "cpptraj_cacovar"

    # Check for a template PDB, if none copy it over
    qa.manage.check_file("template.pdb", "../template.pdb")

    # Check if a pdb trajectory was generated
    file_name = "all_coors.pdb"
    isExist = os.path.isfile(file_name)
    if not isExist:
        print(f"> Can't find the {file_name}.")
        print(f"> Generating it...")
        qa.process.xyz2pdb_traj()

    # Move to the analysis folder
    primary = os.getcwd()
    qa.manage.check_folder("Analysis")
    os.chdir("Analysis")
    qa.manage.check_folder(f"1_{job}")
    os.chdir(f"1_{job}")

    # Copy stored scripts
    qa.manage.copy_script(f"{job}.in")
    qa.manage.copy_script(f"{job}.sh")

    # Execute the script if does not exist or the settings request a recompute
    file_name = "cacovar.dat"
    isExist = os.path.isfile(file_name)
    if not isExist:
        print(f"> Can't find the {file_name}.")
        os.system(f"sbatch {job}.sh")
    elif recompute:
        os.system(f"sbatch {job}.sh")

    # Plot the results if the exist
    if isExist:
        pdb_path = f"{primary}/template.pdb"
        residues = qa.process.get_protein_sequence(pdb_path)
        qa.plot.heatmap(data="cacovar.dat", residues=residues, delete=delete, out_file="matrix_geom.png")


def charge_matrix_analysis(delete, recompute=False) -> None:
    """
    Analyze the charging coupling across all amino acids.

    Performs the analyze to generate matrix files.
    Automates plotting of the results.

    Parameters
    ----------
    delete: List[List[int]]
        A list containing two lists. The first is the residues to delete.
        The second is the columns of the matrix to delete.
    recompute: bool
        Recompute the calculation even if results are already present.
    
    """

    # Check for a template PDB, if none copy it over
    qa.manage.check_file("template.pdb", "../template.pdb")

    # Check for required files
    required_file = "all_charges.xls"
    if not os.path.exists(required_file):
        print(f"{required_file} not found. See qa.process.combine_restarts()") 

    # Check if Analysis folder exits, if it doesn't create one
    primary = os.getcwd() # Save this as our primary home directory for later
    analysis_dir = "Analysis"
    qa.manage.check_folder(analysis_dir)
    os.chdir(analysis_dir)
    job_dir = "2_charge_matrix"
    qa.manage.check_folder(job_dir)
    os.chdir(primary)

    # Generate the charge matrix and mutual information data
    prefix = f"{analysis_dir}/{job_dir}/"
    out_files = ["mimatbb.csv", "chargematbb.csv", "chargematsc.csv"]
    if not os.path.exists(f"{prefix}{out_files[0]}") or not os.path.exists(f"{prefix}{out_files[1]}") or not os.path.exists(f"{prefix}{out_files[2]}"):
        print("> Computing the charge matrix. Please wait...")
        charge_matrix()
        out_files = ["mimatbb.csv", "chargematbb.csv", "chargematsc.csv"]
        shutil.move(out_files[0], f"{analysis_dir}/{job_dir}/{out_files[0]}")
        shutil.move(out_files[1], f"{analysis_dir}/{job_dir}/{out_files[1]}")
        shutil.move(out_files[2], f"{analysis_dir}/{job_dir}/{out_files[2]}")
    elif recompute:
        print("> Computing the charge matrix. Please wait...")
        charge_matrix()
        out_files = ["mimatbb.csv", "chargematbb.csv", "chargematsc.csv"]
        shutil.move(out_files[0], f"{analysis_dir}/{job_dir}/{out_files[0]}")
        shutil.move(out_files[1], f"{analysis_dir}/{job_dir}/{out_files[1]}")
        shutil.move(out_files[2], f"{analysis_dir}/{job_dir}/{out_files[2]}")

    # Move working directory to job directory in Analysis for plotting
    os.chdir(f"{analysis_dir}/{job_dir}")

    # plot the results
    pdb_path = f"{primary}/template.pdb"
    residues = qa.process.get_protein_sequence(pdb_path)
    plot_name = "matrix_charge.png"
    data_name = "chargematbb.csv"
    if not os.path.exists(plot_name):
        qa.plot.heatmap(data=data_name, residues=residues, delete=delete, out_file=plot_name)
    elif recompute:
        qa.plot.heatmap(data=data_name, residues=residues, delete=delete, out_file=plot_name)

    # plot the mututal informatino results
    plot_name = "mi_charge.png"
    data_name = "mimatbb.csv"
    if not os.path.exists(plot_name):
        qa.plot.heatmap(data=data_name, residues=residues, delete=delete, out_file=plot_name, cmap="Blues", v=[0,.2])
    elif recompute:
        qa.plot.heatmap(data=data_name, residues=residues, delete=delete, out_file=plot_name, cmap="Blues", v=[0,.2])


def multiwfn_esp():
    """
    Calculates the electrostatic potential of the system using Multiwfn.

    """

    # Charge paradigms to compute with their Multiwfn codes
    charge_schemes = {'Hirshfeld': '1', 'Voronoi':'2', 'Mulliken': '5', 'ADCH': '11'}



if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    charge_matrix_analysis([[],[]], recompute=False)
