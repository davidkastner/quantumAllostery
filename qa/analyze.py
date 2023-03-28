"""Functions for analyzing an AIMD trajectory."""

import os
import sys
import glob
import time
import shutil
import subprocess
import prody
import numpy as np
import pandas as pd
from joblib import parallel_backend
from sklearn.feature_selection import mutual_info_regression
import qa.plot
import qa.manage
import qa.process
import select
from typing import List


def charge_matrix(pdbfile) -> None:
    """
    Generates mutual information and cross-correlation matrices.

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


    # Write out matrixes to csvs 
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
    joint_df: pd.DataFrame()
        A dataframe with the charge information for both correlated residues.

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


def compute_rmsd(matrix_A, matrix_B):
    """
    Computes the RMSD for two sets of atoms.

    Given two pairs of comparable atom sets, compute the RMSD between them.
    Before computing the RMSD we first use prody to align the structures.

    Parameters
    ----------
    matrix_A: np.array()
        A matrix of atoms with columns as the x, y, and z coordinates.
        The rows, are the different atoms.
    matrix_B: np.array()
        A matrix of atoms with columns as the x, y, and z coordinates.
        The rows, are the different atoms.

    Returns
    -------
    rmsd: float
        The computed RMSD in angstoms.

    Notes
    -----
    Prody is the same underlying alignment algorithm as PyMol and Chimera.

    """
    # Check if the two matrices are the same shape and throw an error otherwise
    if matrix_A.shape != matrix_B.shape:
        raise ValueError("Matrices A and B must have the same shape")

    # Calculate the optimal rotation matrix using ProDy's superpose() function
    _, transformation_matrix = prody.calcTransformation(matrix_A, matrix_B)

    # Apply the transformation matrix to matrix_B to align it with matrix_A
    matrix_B_aligned = transformation_matrix.apply(matrix_B)

    # Calculate the RMSD between the aligned structures
    squared_diff = np.square(matrix_A - matrix_B_aligned)
    sum_squared_diff = np.sum(squared_diff, axis=1)
    mean_squared_diff = np.mean(sum_squared_diff)
    rmsd = np.sqrt(mean_squared_diff)

    return rmsd


def get_rmsd(ref_atoms: List[int], traj_atoms: List[int]):
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
                    [list(map(float, frame_lines[i - 1].split()[1:4])) for i in traj_atoms],
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


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    calculate_charge_schemes()
