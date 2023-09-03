"""Functions for analyzing an AIMD trajectory."""

import os
import sys
import glob
import time
import shutil
import subprocess
import numpy as np
import pandas as pd
from joblib import parallel_backend
from sklearn.feature_selection import mutual_info_regression
import qa.plot
import qa.manage
import qa.process
import select
from typing import List
import MDAnalysis as mda
from scipy.spatial.transform import Rotation as R
from scipy.spatial import distance
import io
from MDAnalysis.analysis import distances
from itertools import combinations


def charge_matrix(pdbfile) -> None:
    """
    Generates mutual information and cross-correlation matrices.

    Parameters
    ----------
    pdbfile : str
        The name of the pdbfile

    """

    start_time = time.time()  # Used to report the executation speed
    # Set variables
    old_res_num = 0
    reslist = []
    atw_all_list = []
    atwsclist = []

    # Search for the reference PDB
    with open(pdbfile, "r") as pdbfile:
        pdbfile_lines = pdbfile.readlines()

        for line in pdbfile_lines:
            if "ATOM" not in line:
                continue
            line_pieces = line.split()
            atom_num = line_pieces[1]
            atom = line_pieces[2]
            res_name = line_pieces[3]
            res_num = int(line_pieces[4])

            # Start new residue
            if res_num != old_res_num:
                reslist.append(res_name)
                old_res_num = res_num
                atw_all_list.append([atom_num])
                atwsclist.append([])
            else:
                atw_all_list[old_res_num - 1].append(atom_num)

            # Only side chain atoms in the atwsclist
            if atom not in ["N", "O", "C", "H"]:
                atwsclist[old_res_num - 1].append(atom_num)

    with open("all_charges.xls", "r") as charge_file:
        charge_file_lines = charge_file.readlines()
        all_chargearray = []
        scchargearray = []
        charge_file_length = len(charge_file_lines)
        atw_all_list_len = len(atw_all_list)
        print(f"> {round(charge_file_length * 0.0005, 2)} ps collected so far")

        for line2 in range(1, charge_file_length):
            atom_charges = charge_file_lines[line2].split()
            all_chargearray.append([])
            scchargearray.append([])

            for residue_index in range(0, atw_all_list_len):
                rescharge = 0.00
                # sum up all the charges for the current residue
                for atom_number in atw_all_list[residue_index]:
                    rescharge += float(
                        atom_charges[int(atom_number) - 1]
                    )
                all_chargearray[line2 - 1].append(rescharge)

                ressccharge = 0.00
                # sum up all the charges for the current residue
                for atom_number in atwsclist[residue_index]:
                    ressccharge += float(
                        atom_charges[int(atom_number) - 1]
                    )
                scchargearray[line2 - 1].append(ressccharge)

    print("> Finished looping through charges")

    reslist_len = len(reslist)
    all_chargearray = np.array(all_chargearray)
    scchargearray = np.array(scchargearray)
    all_chargemutinf = all_chargearray.transpose()
    MI_mat = []

    print(f"> Start looping through {reslist_len} residues")

    for i in range(0, reslist_len):
        print(f"   > rowMI {i}")
        # mutual_info_regression with the transpose matrix of 3k+ charges to 30 residues
        # and the associated row of residues
        rowMI = mutual_info_regression(all_chargearray, all_chargemutinf[i])
        MI_mat.append(rowMI)

    corrcoef = np.corrcoef(all_chargemutinf)
    corrcoefsc = np.corrcoef(scchargearray.transpose())

    # Construct strings from matrices for writing to csv
    mimat_string = ""
    chargemat_all_string = ""
    chargematsc_string = ""
    for resx in range(0, reslist_len):
        for resy in range(0, reslist_len):
            extra = "" if resy == reslist_len - 1 else ","
            mimat_string += str(MI_mat[resx][resy]) + extra
            chargemat_all_string += str(corrcoef[resx][resy]) + extra
            chargematsc_string += str(corrcoefsc[resx][resy]) + extra

        mimat_string += "\n"
        chargemat_all_string += "\n"
        chargematsc_string += "\n"


    # Write out matrixes to CSV's 
    with open("mimatbb.csv", "w") as mimat:
        mimat.write(mimat_string)

    with open("chargematbb.csv", "w") as chargemat:
        chargemat.write(chargemat_all_string)

    with open("chargematsc.csv", "w") as chargemat:
            chargemat.write(chargematsc_string)


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


def get_joint_qres(res_x, res_y, axes_range):
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
    joint_df: pd.DataFrame()
        A dataframe with the charge information for both correlated residues.

    Notes
    -----
    Big picture flow of function.
        1. Read charges.xls in as a pd dataframe
        2. Get atom indices for requested atoms from get_res_atom_indices()
        3. Extract and sum residue columns
        4. Save as a csv
        5. Return as a pandas dataframe
        6. Use joint_plot() to plot the results

    """
    start_time = time.time()  # Used to report the executation speed
    # Create folder for the analysis
    root = os.getcwd()
    os.chdir("Analysis")
    if not os.path.isdir("3_coupling"):
        os.makedirs("3_coupling")
    os.chdir(root)

    #
    structure = 0

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

        # Add the summed residue to the new dataframe
        joint_df[res] = summed_charges

    # Save out the dataframe
    joint_df.to_csv("charge_df.csv")

    extensions = ["png", "svg"]
    for ext in extensions:
        plot_name = f"{res_x}_{res_y}_dist.{ext}"
        qa.plot.get_charge_distributions(joint_df, plot_name, res_x, res_y, ext, axes_range)

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
        qa.plot.heatmap(
            data="cacovar.dat",
            residues=residues,
            delete=delete,
            out_file="matrix_geom.png",
        )


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
    primary = os.getcwd()  # Save this as our primary home directory for later
    analysis_dir = "Analysis"
    qa.manage.check_folder(analysis_dir)
    os.chdir(analysis_dir)
    job_dir = "2_charge_matrix"
    qa.manage.check_folder(job_dir)
    os.chdir(primary)

    # Generate the charge matrix and mutual information data
    prefix = f"{analysis_dir}/{job_dir}/"
    out_files = ["mimatbb.csv", "chargematbb.csv", "chargematsc.csv"]
    if (
        not os.path.exists(f"{prefix}{out_files[0]}")
        or not os.path.exists(f"{prefix}{out_files[1]}")
        or not os.path.exists(f"{prefix}{out_files[2]}")
    ):
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
        qa.plot.heatmap(
            data=data_name, residues=residues, delete=delete, out_file=plot_name
        )
    elif recompute:
        qa.plot.heatmap(
            data=data_name, residues=residues, delete=delete, out_file=plot_name
        )

    # plot the mututal informatino results
    plot_name = "mi_charge.png"
    data_name = "mimatbb.csv"
    if not os.path.exists(plot_name):
        qa.plot.heatmap(
            data=data_name,
            residues=residues,
            delete=delete,
            out_file=plot_name,
            cmap="Blues",
            v=[0, 0.2],
        )
    elif recompute:
        qa.plot.heatmap(
            data=data_name,
            residues=residues,
            delete=delete,
            out_file=plot_name,
            cmap="Blues",
            v=[0, 0.2],
        )


def calculate_charge_schemes():
    """
    Calculated a variety of charge schemes with Multiwfn.

    Compute metal-centered Hirshfeld, Voronoi, Mulliken, and ADCH charges.
    Uses the Multiwfn package to compute the charge partitions.

    Parameters
    ----------
    molden: str
        The name of the molden file to be processed with Multiwfn

    """
    # Get the prefix of the molden file
    molden_files = glob.glob("*.molden")
    if len(molden_files) == 1:
        molden = molden_files[0]
        molden = molden.split(".")[0]
    elif len(molden_files) == 0:
        raise Exception("No molden was found.")
    else:
        raise Exception("More than one molden was found.")

    # If installed correctly, Multiwfn can be called with Multiwfn
    threads = 4
    # Setting the threads too high can cause peformance problems
    command = f"Multiwfn {molden}.molden -nt {threads}"
    print(f"> Using {threads} threads for Multiwfn")

    start_time = time.time()  # Used to report the executation speed
    job_count = 0  # Keep track of the number of calculations to inform user
    # Different charge schemes to be calculated with Multiwfn
    charge_schemes = {"Hirshfeld": "1", "Voronoi": "2", "Mulliken": "5", "ADCH": "11"}
    # Store extracted multiwfn parameters
    temp_dict = {}

    for scheme in charge_schemes:
        print(f"> Current charge scheme: {scheme}")
        proc = subprocess.Popen(
            command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True
        )
        calc_command = charge_schemes[scheme]

        # For atomic charge type corresponding to dict key
        commands = ["7", calc_command, "1", "y", "0", "0", "q"]

        # Run Multiwfn
        output = proc.communicate("\n".join(commands).encode())

        # Watch the Multiwfn output in real time for troubleshooting
        # proc.stdin.write("\n".join(commands).encode())
        # proc.stdin.close()

        # verbose = True
        # if verbose:
        #     for line in iter(proc.stdout.readline, b''):
        #         print(f"   > {line.decode().strip()}")

        # Process the output of Multiwfn
        lines = str(output[0]).split("\\n")
        start = False
        counter = 0

        for num, line in enumerate(lines):
            print(line)
            if "Total valences and free valences" in line:
                start = True
                continue
            if start and "Fe" in line:
                temp_dict[d] = {
                    "Total valence": float(line.split()[-2]),
                    "Free valence": float(line.split()[-1]),
                }
                counter += 1

        # Rename the file as the new, desired name
        new_name = f"{molden}_{scheme}.txt"
        os.rename(f"{molden}.chg", new_name)

    job_count += 1

    total_time = round(time.time() - start_time, 3)  # Time to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Calculated charge schemes for {job_count} QM jobs.
        \tOUTPUT: Generated charge schemes in the current directory.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


def calculate_esp(component_atoms, scheme):
    """
    Calculate the electrostatic potential (ESP) of a molecular component.

    Takes the output from a Multiwfn charge calculation and calculates the ESP.
    Run it from the folder that contains all replicates.
    It will generate a single csv file with all the charges for your residue,
    with one component/column specified in the input residue dictionary.

    Parameters
    ----------
    component_atoms: List[int]
        A list of the atoms in a given component

    """
    # Physical constants
    k = 8.987551 * (10**9)  # Coulombic constant in kg*m**3/(s**4*A**2)
    A_to_m = 10 ** (-10)
    KJ_J = 10**-3
    faraday = 23.06  # Kcal/(mol*V)
    C_e = 1.6023 * (10**-19)
    one_mol = 6.02 * (10**23)
    cal_J = 4.184
    component_esp_list = []

    # Open a charge scheme file as a pandas dataframe
    file_path = glob.glob(f"*_{scheme}.txt")[0]
    df_all = pd.read_csv(file_path, sep="\s+", names=["Atom", "x", "y", "z", "charge"])

    # The index of the metal center assuming iron Fe
    metal_index = df_all.index[df_all["Atom"] == "Fe"][0]
    component_atoms.append(metal_index)

    # Select rows corresponding to an atoms in the component
    df = df_all[df_all.index.isin(component_atoms)]
    df.reset_index(drop=True, inplace=True)

    # Get the new index of the metal as it will have changed
    metal_index = df.index[df["Atom"] == "Fe"][0]

    # Convert columns lists for indexing
    atoms = df["Atom"]  # Now contains only atoms in component
    charges = df["charge"]
    xs = df["x"]
    ys = df["y"]
    zs = df["z"]

    # Determine position and charge of the target atom
    xo = xs[metal_index]
    yo = ys[metal_index]
    zo = zs[metal_index]
    chargeo = charges[metal_index]
    total_esp = 0

    for idx in range(0, len(atoms)):
        if idx == metal_index:
            continue
        else:
            # Calculate esp and convert to units (A to m)
            r = (
                ((xs[idx] - xo) * A_to_m) ** 2
                + ((ys[idx] - yo) * A_to_m) ** 2
                + ((zs[idx] - zo) * A_to_m) ** 2
            ) ** 0.5
            total_esp = total_esp + (charges[idx] / r)

    # Note that cal/kcal * kJ/J gives 1
    component_esp = k * total_esp * ((C_e)) * cal_J * faraday

    return component_esp


def compute_rmsd(matrix_A: np.ndarray, matrix_B: np.ndarray) -> float:
    """
    Computes the Root Mean Square Deviation (RMSD) between two matrices of atoms with their xyz coordinates.

    Parameters
    ----------
    matrix_A : np.ndarray
        A matrix of size N x 3, where N is the number of atoms and the columns are the x, y, and z coordinates.
    matrix_B : np.ndarray
        A matrix of size N x 3, where N is the number of atoms and the columns are the x, y, and z coordinates.

    Returns
    -------
    rmsd : float
        The computed RMSD between the two matrices.
    """
    # Compute the centroids of each matrix
    centroid_A = np.mean(matrix_A, axis=0)
    centroid_B = np.mean(matrix_B, axis=0)

    # Center each matrix by subtracting their respective centroids
    centered_A = matrix_A - centroid_A
    centered_B = matrix_B - centroid_B

    # Compute the covariance matrix
    covariance_matrix = np.dot(centered_A.T, centered_B)

    # Compute the optimal rotation matrix using Singular Value Decomposition (SVD)
    U, _, Vt = np.linalg.svd(covariance_matrix)
    optimal_rotation_matrix = np.dot(U, Vt)

    # Check for reflections
    if np.linalg.det(optimal_rotation_matrix) < 0:
        Vt[-1, :] *= -1
        optimal_rotation_matrix = np.dot(U, Vt)

    # Rotate matrix B using the optimal rotation matrix
    rotated_B = np.dot(centered_B, optimal_rotation_matrix.T)

    # Compute the RMSD
    rmsd = np.sqrt(np.mean(np.sum((centered_A - rotated_B)**2, axis=1)))

    return rmsd


def get_rmsd(ref_atoms: List[int], traj_atoms: List[List[int]]):
    """
    Gets RMSD for specific atoms for different analogs.

    Loops over analog directories and computes a series of RMSDs,
    between a reference set of atoms and a set of atoms from an xyz trajectory.
    The atom indices for the reference structure are given with ref_atoms.
    The atom's xyz coordinates for a frame in an xyz trajectory are traj_atoms.
    ref_atoms and traj_atoms are used to generate matrices A and B.
    The RMSD is then computed with compute_rmsd.

    Parameters
    ----------
    ref_atoms: List[int]
        A list of the atom indices corresponding to the reference xyz structure
    traj_atoms: List[int]
        A list of the atom indices corresponding to the current xyz frame,
        which we will compare the reference atoms to

    Returns
    -------
    rmsd_list: List[List[float]]
        List of lists where each list represents and analog,
        and each list contains RMSDs for each frame.

    Notes
    -----
    Run from the folder with multiple analogs, e.g., MC6, MC6*, MC6*a

    See Also
    --------
    qa.analyze.compute_rmsd()

    """
    start_time = time.time()  # Used to report the executation speed
    frame_count = 0 # Used to report back to the user

    # Input files
    traj_xyz = "all_coors.xyz"
    ideal_xyz = "reference.xyz"

    # Store RMSDs for each analog in a new list
    rmsd_list: List[List[int]] = []

    # Get the directories of each replicate
    primary = os.getcwd()
    analogs = sorted(glob.glob("*/"))

    # Extract reference atoms from ideal_xyz and create matrix_A
    with open(ideal_xyz, "r") as ideal:
        ideal_lines = ideal.readlines()[2:]  # Skip first two lines (number of atoms and comment)
        matrix_A = np.array(
            [list(map(float, ideal_lines[i - 1].split()[1:4])) for i in ref_atoms],
            dtype=float
        )

    # Loop through the analogs
    for index, analog in enumerate(analogs):
        os.chdir(analog)

        # Extract trajectory atoms from traj_xyz and create matrix_B for each frame
        frame_rmsds = []
        with open(traj_xyz, "r") as traj:
            traj_lines = traj.readlines()
            num_atoms = int(traj_lines[0])  # Read the number of atoms from the first line
            num_frames = len(traj_lines) // (num_atoms + 2)  # Calculate the number of frames

            # Loop over frames
            for frame in range(num_frames):
                offset = 2 + frame * (num_atoms + 2)  # Calculate the starting line number for the current frame
                frame_lines = traj_lines[offset:offset + num_atoms]  # Extract the current frame
                matrix_B = np.array(
                    [list(map(float, frame_lines[i - 1].split()[1:4])) for i in traj_atoms[index]],
                    dtype=float
                )
                
                # Compute the RMSD for each frame
                rmsd = compute_rmsd(matrix_A, matrix_B)
                frame_rmsds.append(rmsd)
                frame_count += 1

        rmsd_list.append(frame_rmsds)
        os.chdir(primary)

    total_time = round(time.time() - start_time, 3)
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tOUTPUT: Computed the RMSD for {index + 1} analogs.
        \tOUTPUT: Computed the RMSD for {frame_count} frames.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )

    return rmsd_list


def td_coupling(res_x, res_y, replicate_dir):
    """
    Calculates and plots the partial charge for two residues over time.

    This function will calculate the charge for each frame for two residues,
    and plot them. This will help show if the charges are inversely correlated.

    Parameters
    ----------
    res_x : str
        The first amino acid to be plotted against time.
    res_y : str
        The second amino acid to be plotted against time.
    replicate_dir: str
        The name of the replicate that we will calculate the time-dependent
        coupling for

    Returns
    -------
    joint_df: pd.DataFrame()
        A dataframe with the charge information for both correlated residues.

    See Also
    --------
    qa.analyze.get_joint_qres()

    """
    start_time = time.time()  # Used to report the executation speed
    # Create folder for the results in Analysis
    root = os.getcwd()
    os.chdir(f"{replicate_dir}/Analysis")
    if not os.path.isdir("4_time_coupling"):
        os.makedirs("4_time_coupling")
    os.chdir(root)

    # Create a new DataFrame to append the two residues of interest
    residues = [res_x, res_y]
    joint_df = pd.DataFrame(columns=residues)

    # Get the indices of the atoms to parse the charge.xls file
    os.chdir(replicate_dir)
    for res in residues:
        # Open the combined charge file and read it in
        charge_file = "all_charges.xls"
        charge_df = pd.read_csv(charge_file, sep="\t")

        # Get the atom indices of the residue
        atom_indices = qa.process.get_res_atom_indices(res)
        # Sum the charges of the atoms of the requested residues
        summed_charges = charge_df[charge_df.columns[atom_indices]].sum(axis=1).tolist()

        # Add the summed residue to the new dataframe
        joint_df[res] = summed_charges

    ext = "svg"
    plot_name = f"{res_x}_{res_y}_td.{ext}"
    qa.plot.time_coupling_plot(joint_df, plot_name, res_x, res_y, ext)

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t-------------------------GET TD CHARGES END-------------------------
        \tRESULT: Extracted and computed the joint charge distribution.
        \tOUTPUT: Created the charge distribution plot: {plot_name}.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )

    return joint_df

def pairwise_distances_csv(pdb_traj_path, output_file):
    """
    Calculate pairwise distances between residue centers of mass and save the result to a CSV file.
    
    Parameters
    ----------
    pdb_traj_path : str
        The file path of the PDB trajectory file.
    output_file : str
        The name of the output CSV file.
    """

    # Read the trajectory file and split it into models
    with open(pdb_traj_path) as f:
        models = f.read().split("END")
    
    # Create a list of StringIO objects for each model
    frame_files = [io.StringIO(model) for model in models if model.strip()]
    universes = [mda.Universe(frame_file, format="pdb") for frame_file in frame_files]

    # Generate column names based on residue pairs
    residue_names = [residue.resname + str(residue.resid) for residue in universes[0].residues]
    residue_pairs = list(combinations(residue_names, 2))
    column_names = [f"{pair[0]}-{pair[1]}" for pair in residue_pairs]

    pairwise_distances = []
    for universe in universes:
        # Calculate the center of mass for each residue
        residue_com = np.array([residue.atoms.center_of_mass() for residue in universe.residues])
        
        # Calculate the pairwise distance matrix
        distance_matrix = distances.distance_array(residue_com, residue_com)
        pairwise_distances.append(distance_matrix[np.triu_indices(len(residue_com), k=1)])

    # Create a DataFrame with pairwise distances and column names
    pairwise_distances_df = pd.DataFrame(pairwise_distances, columns=column_names)
    pairwise_distances_df.to_csv(output_file, index=False)

def parse_components(components):
    """
    Parse the components input, handling range inputs (e.g., '253-424').

    Parameters
    ----------
    components : list of lists
        A list of components, each component being a list of strings representing atom numbers or ranges

    Returns
    -------
    parsed_components : list of lists
        A list of components, each component being a list of integers representing atom numbers
    """
    parsed_components = []
    for component in components:
        parsed_component = []
        for atom in component:
            if '-' in atom:  # Handle range input
                start, end = map(int, atom.split('-'))
                parsed_component.extend(list(range(start, end + 1)))
            else:
                parsed_component.append(int(atom))
        parsed_components.append(parsed_component)
    return parsed_components

def centroid_distance(components):
    """
    Calculate the distance between two structural components.

    This purpose of this script is to do more than calculate a simple distance.
    It will take in two sets of amino acids.
    This can be a single amino acid or an arbitrary number of amino acid.
    It will then calculate the centroid of the two components,
    and the distance between the centroids.

    Parameters
    ----------
    components : List of lists
        A list of two lists corresponding to two components with their atoms e.g., [[487],[253-424]]

    Notes
    -----
    Examples of interesting components for mimochromes:

        all : 1-487
        lower : 1-252
        upper : 253-424
        lower-his : 1-86,104-252
        heme : 425-486
        his : 87-103

    """
    # Handle range inputs in components
    components = parse_components(components)

    start_time = time.time()  # Used to report the executation speed
    ignore = ["Analysis/", "coordinates/", "inputfiles/"]
    structure_count = 0

    # Directory containing all replicates
    primary_dir = os.getcwd()
    directories = sorted(glob.glob("*/"))
    replicates = [i for i in directories if i not in ignore]
    replicate_count = len(replicates)  # Report to user

    # Create a pandas dataframe to store the centroids and the distances
    centroid_dist_df = pd.DataFrame(columns=["centroid_1x", "centroid_1y", "centroid_1z", 
                                             "centroid_2x", "centroid_2y", "centroid_2z", 
                                             "distance"])

    for replicate in replicates:
        os.chdir(replicate)
        secondary_dir = os.getcwd()

        # Change into the directory with the structures
        os.chdir("coordinates")
        structures = sorted(glob.glob("*.xyz"), key=lambda x: int(os.path.splitext(x)[0]))

        for structure in structures:
            with open(structure, "r") as structure_coords:
                structure_lines = structure_coords.readlines()[2:]  # Skip the first two lines

                centroids = []
                for component in components:
                    # Parse the XYZ coordinates for the selected atoms
                    selected_atom_coords = [structure_lines[i-1].split()[1:] for i in component]
                    selected_atom_coords = np.array(selected_atom_coords, dtype=float)

                    # Calculate the centroid of each component
                    centroid = np.mean(selected_atom_coords, axis=0)
                    centroids.append(centroid)

                # Calculate the distance between the two centroids
                dist = distance.euclidean(centroids[0], centroids[1])

                # Save the centroids and the distance to centroid_dist_df
                new_row = pd.DataFrame({
                    "centroid_1x": [centroids[0][0]], "centroid_1y": [centroids[0][1]], "centroid_1z": [centroids[0][2]],
                    "centroid_2x": [centroids[1][0]], "centroid_2y": [centroids[1][1]], "centroid_2z": [centroids[1][2]],
                    "distance": [dist]
                })

                centroid_dist_df = pd.concat([centroid_dist_df, new_row], ignore_index=True)
                structure_count += 1

        # Move back to the replicate directory
        os.chdir(primary_dir)

    # Save the dataframe to a csv file
    centroid_dist_df.to_csv("centroid_distance.csv", index=False)

    total_time = round(time.time() - start_time, 3)  # Time to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Performed operation on {replicate_count} replicates.
        \tOUTPUT: Computed {structure_count} centroid distances.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )



if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    calculate_charge_schemes()
