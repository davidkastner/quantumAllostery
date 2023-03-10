"""Process the raw TeraChem output for future analysis."""

import os
import re
import sys
import glob
import time
import shutil
import subprocess
import numpy as np
import pandas as pd
from collections import OrderedDict
from typing import List, Tuple
from biopandas.pdb import PandasPdb
import qa.reference


def get_pdb() -> str:
    """
    Searches all directories recursively for a PDB file.

    If more than one PDB is found it will use the first one.
    If no PDB file was found, it will prompt the user for a PDB file path.

    Returns
    -------
    pdb_file : str
        The path of a PDB file within the current directory (recursive).

    Notes
    -----
    Currently it uses the name to distinguish single structures,
    ensembles, and trajectories.
    In the future, this function should check the contents to confirm.

    """
    # A list of all PDB's found after recursive search
    pdb_files = sorted(glob.glob("./**/template.pdb", recursive=True))
    for index, pdb in enumerate(pdb_files):
        # Trajectory PDB's should be marked as ensemble or traj
        if "ensemble" in pdb or "traj" in pdb or "top" in pdb:
            continue
        else:
            pdb_file = pdb_files[index]

    # Multiple or no PDB files found scenarios
    if len(pdb_files) == 1:
        print(f"> Using {pdb_file} as the template PDB.")
    elif len(pdb_files) > 1:
        print(f"> More than one PDB file found -> Using {pdb_file}.")
    else:
        pdb_file = input("No PDB files was found. What is the path to your PDB file? ")

    return pdb_file


def get_xyz() -> str:
    """
    Searches all directories for a XYZ file.

    If more than one XYZ is found it will use the first one.
    If no XYZ file was found it will prompt the user for the XYZ file path.

    Returns
    -------
    xyz_name : str
        The path of a XYZ file within the current directory.

    """
    # Search recursively for an xyz file
    xyz_names = glob.glob("./**/*.xyz", recursive=True)

    # Check the results to confirm that there was only one xyz file found
    if len(xyz_names) == 1:
        xyz_name = xyz_names[0]
        print(f"> Found the coordinate file {xyz_name}.")
    elif len(xyz_names) > 1:
        xyz_name = xyz_names[0]
        print(f"> Found more than one XYZ, using {xyz_name}.")
    else:
        xyz_name = input("No XYZ was found. What is the path to your XYZ? ")

    return xyz_name


def get_atom_count() -> int:
    """
    Finds an xyz file and gets the number of atoms.

    Returns
    -------
    atom_count : int
        The number of atoms in the identified xyz file.
    """
    # Find an xyz file
    xyz_name = get_xyz()
    with open(xyz_name, "r") as xyz_file:
        atom_count = int(xyz_file.readline().strip())

    return atom_count


def get_protein_sequence(pdb_path) -> List[str]:
    """
    Gets the full amino acid sequence of your protein.

    See Also
    --------
    qa.plot.heatmap()

    """

    # Get the template file and load it as a pandas dataframe
    pdb_df = PandasPdb().read_pdb(pdb_path).df["ATOM"]

    # Filter the dataframe so there is one entry for each residue
    residues_df = pdb_df[["residue_name", "residue_number"]]
    residues_df = residues_df.drop_duplicates(subset=["residue_number"], keep="first")

    # Convert it to a list of amino acids
    residues_df = residues_df["residue_name"]
    residues_list = residues_df.values.tolist()

    return residues_list


def get_charge_file() -> str:
    """
    Searches all directories for a charge xls file.

    If more than one .xls file is found it will use the first one.
    If no .xls file was found it will prompt the user for the .xls file path.
    This is the standard charge output for TeraChem.

    Returns
    -------
    charge_file : str
        The path of a charge .xls file within the current directory.

    """

    # Get the xls from the current directory
    charge_files = sorted(glob.glob("*.xls"))

    # Check the results to confirm that there was only one .xls file found
    if len(charge_files) == 1:
        charge_file = charge_files[0]
        print(f"> Found the charge file {charge_file}.")
    elif len(charge_files) > 1:
        charge_file = charge_files[0]
        print(f"> Found more than one .xls, using {charge_file}.")
    else:
        charge_file = input("> No .xls was found. What is the path to your .xls? ")

    return charge_file


def combine_restarts(
    atom_count, all_charges: str = "all_charges.xls", all_coors: str = "all_coors.xyz"
) -> None:
    """
    Collects all charges or coordinates into single xls and xyz files.

    Likely the first executed function after generating the raw AIMD data.
    Trajectories were likely generated over multiple runs.
    This function combines all coordinate and charge data for each run.

    Parameters
    ----------
    all_charges : str
        The name of the file containing all charges in xls format.
    all_coors.xyz : str
        The name of the file containing the coordinates in xyz format.
    atom_count : int
        The number of atoms in the structure

    Notes
    -----
    Run from the directory that contains the run fragments.

    See Also
    --------
    combine_replicates: Combines all the combined trajectories of each replicate into one.

    """

    start_time = time.time()  # Clock executation speed

    # Collect all qmscript.out files
    out_files = glob.glob("./**/qmscript.out", recursive=True)
    out_files_count = len(out_files)  # Will report this to user
    # run_info: The MD restart step, the outfile name, and its scr directory
    run_info: list[list[int, str, str]] = []

    # Get the first MD step from out_files to identify where each job restarted
    # Note: Restarting occurs from the last restart point not the last frame
    for out_file in out_files:
        out_content = open(out_file, "r").readlines()
        for line in out_content:
            if "MD STEP" in line:
                # A list of info about each run [start step, file, scr dir]
                md_step = int(line.split()[4])
                run_info.append([md_step, out_file])
                break

    # Get the src directory locations for the restart attempt to generate paths
    scrdir = glob.glob("./**/scr*/", recursive=True)
    scrdir.sort()  # Sort by modification date
    run_info.sort()  # Sort by MD step
    for index, step in enumerate(run_info):
        # Add the name of the scr directory location to the run_info list
        step.append(scrdir[index])

    # Delete so we don't append to a previous version when rerunning
    if os.path.exists(all_charges):
        os.remove(all_charges)
    if os.path.exists(all_coors):
        os.remove(all_coors)

    # Use the run_info to open the charge and coordinate files
    first_run = True
    for index, run in enumerate(run_info):
        coors_file = open(f"{run[2]}coors.xyz", "r").readlines()
        charge_file = open(f"{run[2]}charge.xls", "r").readlines()
        # Create combined charge and coors files
        all_coors_file = open(all_coors, "a")
        all_charges_file = open(all_charges, "a")

        # First run
        if first_run:
            coor_run_end = (run_info[index + 1][0] - 1) * (atom_count + 2)
            all_coors_file.writelines(coors_file[:coor_run_end])
            # The first charges line is a special header line that we add once
            # so we don't need to substract one because it cancels with the index offset
            charge_run_end = run_info[index + 1][0]
            all_charges_file.writelines(charge_file[:charge_run_end])
            first_run = False

        # Last run
        elif index == len(run_info) - 1:
            # Go all the way to the end so a frame isn't left off
            all_coors_file.writelines(coors_file)
            all_charges_file.writelines(charge_file[1:])

        # Other run
        else:
            # To get the number of frames in run two,
            # substract the restart number from the frames in run 1
            coor_run_end = ((run_info[index + 1][0] - 1) - (run[0] - 1)) * (
                atom_count + 2
            )
            all_coors_file.writelines(coors_file[:coor_run_end])
            charge_run_end = (run_info[index + 1][0]) - (run[0] - 1)
            all_charges_file.writelines(charge_file[1:charge_run_end])

    # Close files
    all_coors_file.close()
    all_charges_file.close()

    # Check the number of frames in the combined xyz trajectory and print for user
    coors_frame_count = 0
    with open(all_coors, "r") as coors:
        for line in coors:
            if "frame" in line:
                coors_frame_count += 1

    # Check number of charge frames and print for user
    charge_frame_count = -1
    with open(all_charges, "r") as charges:
        for line in charges:
            charge_frame_count += 1

    total_time = round(time.time() - start_time, 3)  # Seconds to run
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: {coors_frame_count} frames and {charge_frame_count} charges from {out_files_count} runs.
        \tOUTPUT: Generated {all_charges} and {all_coors}.
        \tOUTPUT: {os.path.abspath(os.getcwd())}.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


def combine_replicates(
    all_charges: str = "all_charges.xls", all_coors: str = "all_coors.xyz"
) -> None:
    """
    Collects charges or coordinates into a xls and xyz file across replicates.

    Parameters
    ----------
    all_charges : str
        The name of the file containing all charges in xls format.
    all_coors.xyz : str
        The name of the file containing the coordinates in xyz format.

    Notes
    -----
    Run from the directory that contains the replicates.
    Run combine_restarts first for if each replicated was run across multiple runs.
    Generalized to combine any number of replicates.

    See Also
    --------
    combine_restarts: Combines restarts and should be run first.
    """

    # General variables
    start_time = time.time()  # Used to report the executation speed
    files = [all_charges, all_coors]  # Files to be concatonated
    charge_files: list[str] = []  # List of the charge file locations
    coors_files: list[str] = []  # List of the coors file locations
    root = os.getcwd()
    dirs = sorted(glob.glob(f"{root}/*/"))  # glob to efficiently grab only dirs
    replicates = len(dirs)  # Only used to report to user

    # Loop through all directories containing replicates
    for dir in dirs:
        if os.path.isfile(f"{dir}{files[0]}") and os.path.isfile(f"{dir}{files[1]}"):
            charge_files.append(f"{dir}{files[0]}")
            coors_files.append(f"{dir}{files[1]}")

    new_file_names = [f"raw_{all_charges}", all_coors]
    file_locations = [charge_files, coors_files]
    # Loop over the file names and their locations
    for file_name, file_location in zip(new_file_names, file_locations):
        # Open a new file where we will write the concatonated output
        with open(file_name, "wb") as outfile:
            for loc in file_location:
                with open(loc, "rb") as infile:
                    shutil.copyfileobj(infile, outfile)

    # The combined charge file now has multiple header lines
    first_line = True
    with open(new_file_names[0], "r") as raw_charge_file:
        with open(files[0], "w") as clean_charge_file:
            for line in raw_charge_file:
                # We want the first line to have the header
                if first_line == True:
                    clean_charge_file.write(line)
                    first_line = False
                # After the first, no lines should contain atom names
                else:
                    if "H" in line:
                        continue
                    else:
                        clean_charge_file.write(line)
    # Delete the charge file with the extra headers to keep the dir clean
    os.remove(new_file_names[0])

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Combined {replicates} replicates.
        \tOUTPUT: Generated {files[0]} and {files[1]} in the current directory.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


def get_residue_identifiers(template, by_atom=True) -> List[str]:
    """
    Gets the residue identifiers such as Ala1 or Cys24.

    Returns either the residue identifiers for every atom, if by_atom = True
    or for just the unique amino acids if by_atom = False.

    Parameters
    ----------
    template: str
        The name of the template pdb for the protein of interest.
    by_atom: bool
        A boolean value for whether to return the atom identifiers for all atoms

    Returns
    -------
    residues_indentifier: List(str)
        A list of the residue identifiers

    """
    # Get the residue number identifiers (e.g., 1)
    residue_number = (
        PandasPdb().read_pdb(template).df["ATOM"]["residue_number"].tolist()
    )
    # Get the residue number identifiers (e.g., ALA)
    residue_name = PandasPdb().read_pdb(template).df["ATOM"]["residue_name"].tolist()
    # Combine them together
    residues_indentifier = [
        f"{name}{number}" for number, name in zip(residue_number, residue_name)
    ]

    # Return only unique entries if the user sets by_atom = False
    if not by_atom:
        residues_indentifier = list(OrderedDict.fromkeys(residues_indentifier))

    return residues_indentifier


def average_charge_residues(charge_file, template):
    """
    Averages charge file by residue.

    Reduces inaccuracies introduced by the limitations of Mulliken charges.

    Parameters
    ----------
    charge_file: str
        The name of the file containing the charges.
    template: str
        The name of the template pdb for the protein of interest.

    Returns
    -------
    charge_df: pd.DataFrame
        The charge file summed by residue and stored as a pd.DataFrame.

    """
    # Read in the charge data file as pd.DataFrame
    charge_df = pd.read_table(charge_file, sep="\t")
    # charge_df = charge_df.iloc[: , :-1]

    # Get the residue identifiers (e.g., 1Ala) for each atom
    residues_indentifier = get_residue_identifiers(template)

    # Group of the charge data frame by the residue identifiers
    try:
        groups = charge_df.groupby(residues_indentifier, axis=1, sort=False)
    except:
        print("   > ERROR: You likely have an extra hanging column.")
    # Compute the mean of all atoms in a residue
    avg_by_residues = groups.mean()

    return avg_by_residues


def xyz2pdb(xyz_list: List[str]) -> None:
    """
    Converts an xyz file into a pdb file.

    Parameters
    ----------
    xyz_list: List(str)
        A list of file names that you would like to convert to PDB's

    Note
    ----
    Make sure to manually check the PDB that is read in.
    Assumes no header lines.
    Assumes that the only TER flag is at the end.

    """
    start_time = time.time()  # Used to report the executation speed

    # Search for the XYZ and PDB files names
    pdb_name = get_pdb()
    pdb_file = open(pdb_name, "r").readlines()
    max_atom = int(pdb_file[len(pdb_file) - 3].split()[1])

    for index, xyz in enumerate(xyz_list):
        new_file = open(f"{index}.pdb", "w")
        xyz_file = open(xyz, "r").readlines()
        atom = -1  # Start at -1 to skip the XYZ header
        line_count = 0
        for line in xyz_file:
            line_count += 1
            if atom > 0:
                atom += 1
                try:
                    x, y, z = line.strip("\n").split()[1:5]  # Coordinates from xyz file
                except:
                    print(f"> Script died at {line_count} -> '{line}'")
                    quit()
                pdb_line = pdb_file[atom - 2]  # PDB is two behind the xyz
                new_file.write(
                    f"{pdb_line[0:30]}{x[0:6]}  {y[0:6]}  {z[0:6]}  {pdb_line[54:80]}"
                )
            else:
                atom += 1
            if atom > max_atom:
                atom = -1
                new_file.write("END\n")

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Converted {len(xyz_list)} to pdbs.
        \tOUTPUT: Generated in the current directory.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


def xyz2pdb_traj() -> None:
    """
    Converts an xyz trajectory file into a pdb trajectory file.

    Note
    ----
    Make sure to manually check the PDB that is read in.
    Assumes no header lines.
    Assumes that the only TER flag is at the end.

    """

    start_time = time.time()  # Used to report the executation speed

    # Search for the XYZ and PDB files names
    pdb_name = get_pdb()
    xyz_name = get_xyz()

    # Remove the extension to get the protein name to use as the PDB header
    protein_name = pdb_name.split("/")[-1][:-4]
    new_pdb_name = f"all_coors.pdb"

    # Open files for reading
    xyz_file = open(xyz_name, "r").readlines()
    pdb_file = open(pdb_name, "r").readlines()
    max_atom = int(pdb_file[len(pdb_file) - 3].split()[1])
    new_file = open(new_pdb_name, "w")

    atom = -1  # Start at -1 to skip the XYZ header
    line_count = 0
    for line in xyz_file:
        line_count += 1
        if atom > 0:
            atom += 1
            try:
                x, y, z = line.strip("\n").split()[1:5]  # Coordinates from xyz file
            except:
                print(f"> Script died at {line_count} -> '{line}'")
                quit()
            pdb_line = pdb_file[atom - 2]  # PDB is two behind the xyz
            new_file.write(
                f"{pdb_line[0:30]}{x[0:6]}  {y[0:6]}  {z[0:6]}  {pdb_line[54:80]}"
            )
        else:
            atom += 1
        if atom > max_atom:
            atom = -1
            new_file.write("END\n")

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Converted {xyz_name} to {new_pdb_name}.
        \tOUTPUT: Generated {new_pdb_name} in the current directory.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


def xyz2pdb_ensemble() -> None:
    """
    Converts an xyz trajectory file into a pdb ensemble.

    Note
    ----
    Assumes that the only TER flag is at the end.
    """

    # Search for the XYZ and PDB files names
    start_time = time.time()  # Used to report the executation speed
    pdb_name = get_pdb()
    xyz_name = get_xyz()
    # Remove the extension to get the protein name to use as the PDB header
    protein_name = pdb_name.split("/")[-1][:-4].upper()
    new_pdb_name = f"{protein_name}_ensemble.pdb"

    # Open files for reading
    xyz_file = open(xyz_name, "r").readlines()
    pdb_file = open(pdb_name, "r").readlines()
    max_atom = int(pdb_file[len(pdb_file) - 3].split()[1])
    new_file = open(new_pdb_name, "w")

    atom = -1  # Start at -1 to skip the XYZ header
    model_number = 2
    new_file.write(f"{protein_name}\n")  # PDB header line
    new_file.write(f"MODEL        1\n")  # The first line will always be MODEL 1

    for line in xyz_file:
        if atom > 0:
            atom += 1
            x, y, z = line.strip("\n").split()[1:5]  # Coordinates from xyz file
            pdb_line = pdb_file[atom - 2]  # PDB is two behind the xyz
            new_file.write(
                f"{pdb_line[0:30]}{x[0:6]}  {y[0:6]}  {z[0:6]}  {pdb_line[54:80]}\n"
            )
        else:
            atom += 1
        if atom > max_atom:
            atom = -1
            new_file.write("TER\n")
            new_file.write("ENDMDL\n")
            new_file.write(f"MODEL        {str(model_number)}\n")
            model_number += 1

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Converted {xyz_name} to {new_pdb_name}.
        \tOUTPUT: Generated {new_pdb_name} in the current directory.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


def clean_incomplete_xyz() -> None:
    """
    For removing incomplete frames during troublshooting.

    This is current under construction.
    I am not sure about its use cases.

    """

    start_time = time.time()  # Used to report the executation speed
    orig_file = "all_coors.xyz"
    new_file = "all_coors_clean.xyz"
    incomplete = 0  # Only used to create user status report at end

    with open(new_file, "w") as coors_file_new:
        with open(orig_file, "r") as coors_file:
            section_delim = coors_file.readline().strip()  # First line
            section = []  # Stores the lines for each section
            first_line = True  # The first line is a unique case

            for line in coors_file:
                # Write out the first line no matter what
                if first_line:
                    section.append(line)
                    first_line = False
                else:
                    # Reached the end of a section?
                    if line[: len(section_delim)] == section_delim:
                        # Check if the section has all the atoms it should
                        if len(section) == int(section_delim) + 2:
                            # Write the section out to the new file if complete
                            for line in section:
                                coors_file_new.write(line)
                        else:
                            incomplete += 1
                        # Start a new section
                        section = []
                        section.append(line)
                    else:
                        section.append(line)

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Found {incomplete} incomplete sections.
        \tOUTPUT: Generated {new_file} in the current directory.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


def check_valid_resname(res) -> Tuple[str, int]:
    """
    Checks if a valid resname has been identified.

    Excepts a resname of the form e.g. Ala1, A1, Gly12, G12.
    If an incorrect resname is supplied the fuction will exit with an warning.

    Parameters
    ----------
    res : str
        Name of a residue of the form e.g. Ala1, Gly12.

    Return
    ------
    aa_name : str
        The requested amino acid's three letter code.
    aa_name : int
        The requested amino acid's position in the sequence.

    """
    # Get amino acid identifier information from the reference module
    aa_identifiers = qa.reference.get_aa_identifiers()
    # Check if one or three letter code provided by counting letters
    letter_count = sum(map(str.isalpha, res))

    if letter_count == 1:
        aa_name = res[0].upper().strip()  # clean the input
        aa_num = int(res[1:])
        # Check if the provided one-letter code matches a known amino acid
        if not any([True for k, v in aa_identifiers.items() if v[0] == aa_name]):
            raise ValueError("> The provided one-letter code is unknown.")

    if letter_count == 3:
        aa_name = res[:3].upper().strip()  # clean the input
        aa_num = int(res[3:])
        # Check if the provided three-letter code matches a known amino acid
        if not any([True for k, v in aa_identifiers.items() if v[1] == aa_name]):
            raise ValueError("> The provided three-letter code is unknown.")

    print(f"> Reqesting amino acid {aa_name} at index {aa_num}.")

    return aa_name, aa_num


def get_res_atom_indices(res, scheme="all") -> List[int]:
    """
    For a residue get the atom indices of all atoms in the residue.

    Parameters
    ----------
    res : str
        Name of a residue of the form e.g. Ala1, Gly12.
    type :str
        The type of atom indices to retrieve e.g., all, backbone

    Returns
    -------
    residue_indices : list
        A list of all atom indices for a given residue.

    """

    # Function that gets the path of a PDB file
    pdb = get_pdb()
    # Check if the requested resname is valid
    aa_name, aa_num = check_valid_resname(res)
    # Convert the pdb to a pandas dataframe
    ppdb = PandasPdb().read_pdb(pdb).df["ATOM"]

    # Indices for all residues or for just the backbone
    residue_df = ppdb[
        (ppdb["residue_name"] == aa_name) & (ppdb["residue_number"] == aa_num)
    ]
    atom_index_list = residue_df.index.tolist()

    # Use if you only want the backbone atoms summed
    if scheme == "backbone":
        print(
            "> Retrieving only backbone indices. See qa.process.get_res_atom_indices()"
        )
        bb_atoms = ["N", "H", "C", "O"]
        backbone_df = residue_df[residue_df["atom_name"].isin(bb_atoms)]
        atom_index_list = backbone_df.index.tolist()

    if scheme != "all" and scheme != "backbone":
        raise ValueError("> ERROR: Scheme not recognized. Select all or backbone.")

    # Alert the user if the list comes out empty
    if len(atom_index_list) == 0:
        raise ValueError("> ERROR: No atom indices were found. Verify that it exists.")

    return atom_index_list


def clean_qm_jobs(first_job: int, last_job: int, step: int) -> None:
    """
    Cleans all QM jobs and checks for completion.

    We ran single points at a higher level of theory from the SQM simulations.
    Some jobs will inevitable die do to memory or convergence issues.
    It is important to check that all the jobs finished successfully.
    This script checks that all the jobs finished.
    Once it has confirmed that all jobs finished sucessfully,
    it will clean up the QM by deleting log files and scratch directories.

    Parameters
    ----------
    first_job: int
        The name of the first directory and first job e.g., 0
    last_job: int
        The name of the last directory and last job e.g., 39900
    step: int
        The step size between each single point.
    """
    start_time = time.time()  # Used to report the executation speed

    # Directory containing all replicates
    primary_dir = os.getcwd()
    replicates = sorted(glob.glob("*/"))
    total_job_count = 0  # Report to user upon job completion
    incomplete_job_count = 0  # Report to user upon job completion

    for replicate in replicates:
        os.chdir(replicate)
        # The location of the current replicate
        secondary_dir = os.getcwd()
        print(f"> Checking {secondary_dir}.")

        # A list of all job directories assuming they are named as integers
        job_dirs = [str(dir) for dir in range(first_job, last_job, step)]
        for dir in job_dirs:
            total_job_count += 1
            os.chdir(dir)
            tertiary_dir = os.getcwd()

            # Get the out file, we use glob so we can generalize to all out files
            out_name = glob.glob("*.out")
            if len(out_name) < 1:
                print(f"   > Job in {tertiary_dir} did not finish.")
            else:
                # Determine if a job finished
                with open(out_name[0], "r") as out_file:
                    lines = out_file.read().rstrip().splitlines()
                    last_line = lines[-1]
                    # The phrase Job finished is indicative of a success
                    if "Job finished" not in last_line:
                        incomplete_job_count += 1
                        print(f"   > Job in {tertiary_dir} did not finish.")

                    # The job completed, so delete extra scr directories
                    else:
                        scr_dirs = glob.glob("scr*/")
                        # Sort the scr directories by age (oldest to newest)
                        sorted_scr_dirs = sorted(
                            scr_dirs, key=os.path.getmtime, reverse=True
                        )
                        # Only keep the newest
                        for scr_dir in sorted_scr_dirs[1:]:
                            shutil.rmtree(scr_dir)
                            print(f"   > Delete extra scratch directory: {scr_dir}")

            os.chdir(secondary_dir)
        os.chdir(primary_dir)

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Checked {total_job_count} jobs for completion.
        \tOUTPUT: Found {incomplete_job_count} incomplete or problematic jobs.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


def combine_qm_charges(first_job: int, last_job: int, step: int) -> None:
    """
    Combines the charge_mull.xls files generate by TeraChem single points.

    After running periodic single points on the ab-initio MD data,
    we need to process the charge data so that it matches the SQM data.

    Parameters
    ----------
    first_job: int
        The name of the first directory and first job e.g., 0
    last_job: int
        The name of the last directory and last job e.g., 39900
    step: int
        The step size between each single point.

    """
    start_time = time.time()  # Used to report the executation speed
    new_charge_file = "all_charges.xls"
    current_charge_file = "charge_mull.xls"
    ignore = ["Analysis/"]

    # Directory containing all replicates
    primary_dir = os.getcwd()
    directories = sorted(glob.glob("*/"))
    replicates = [i for i in directories if i not in ignore]
    replicate_count = len(replicates)  # Report to user

    for replicate in replicates:
        frames = 1  # Saved to report to the user
        os.chdir(replicate)
        # The location of the current qm job that we are appending
        secondary_dir = os.getcwd()
        print(f"   > Adding { secondary_dir}")

        # Create a new file where we will store the combined charges
        first_charges_file = True  # We need the title line but only once

        if os.path.exists(new_charge_file):
            os.remove(new_charge_file)  # Since appending remove old version
            print(f"      > Deleting old {secondary_dir}/{new_charge_file}.")
        with open(new_charge_file, "a") as combined_charges_file:
            # A list of all job directories assuming they are named as integers
            job_dirs = [str(dir) for dir in range(first_job, last_job, step)]

            # Change into one of the QM job directories
            for index, dir in enumerate(job_dirs):
                os.chdir(dir)
                tertiary_dir = os.getcwd()
                os.chdir("scr")
                # Open an individual charge file from a QM single point
                atom_column = []
                charge_column = []

                # Open one of the QM charge single point files
                with open(current_charge_file, "r") as charges_file:
                    # Separate the atom and charge information
                    for line in charges_file:
                        clean_line = line.strip().split("\t")
                        charge_column.append(clean_line[1])
                        atom_column.append(clean_line[0])

                # Join the data and separate it with tabs
                charge_line = "\t".join(charge_column)

                # For some reason, TeraChem indexes at 0 with SQM,
                # and 1 with QM so we change the index to start at 1
                atoms_line_reindex = []
                for atom in atom_column:
                    atom_list = atom.split()
                    atom_list[0] = str(int(atom_list[0]) - 1)
                    x = " ".join(atom_list)
                    atoms_line_reindex.append(x)
                atom_line = "\t".join(atoms_line_reindex)

                # Append the data to the combined charges data file
                # We only add the header line once
                if first_charges_file:
                    combined_charges_file.write(f"{atom_line}\n")
                    combined_charges_file.write(f"{charge_line}\n")
                    frames += 1
                    first_charges_file = False
                # Skip the header if it has already been added
                else:
                    if "nan" in charge_line:
                        print(f"      > Found nan values in {index * 100}!!")
                    combined_charges_file.write(f"{charge_line}\n")
                    frames += 1

                os.chdir(secondary_dir)
        print(f"      > Combined {frames} frames.")
        os.chdir(primary_dir)

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Combined charges across {replicate_count} replicates.
        \tOUTPUT: Generated {new_charge_file} in the current directory.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )


def combine_qm_replicates() -> None:
    """
    Combines the all_charges.xls files for replicates into a master charge file.

    """
    start_time = time.time()  # Used to report the executation speed
    charge_file = "all_charges.xls"
    ignore = ["Analysis/"]

    # Directory containing all replicates
    primary_dir = os.getcwd()
    directories = sorted(glob.glob("*/"))
    replicates = [i for i in directories if i not in ignore]
    replicate_count = len(replicates)  # Report to user

    # Remove any old version because we are appending
    if os.path.exists(charge_file):
        os.remove(charge_file)
        print(f"      > Deleting old {charge_file}.")

    # Create a new file to save charges
    with open(charge_file, "a") as new_charge_file:
        first_charge_file = True
        for replicate in replicates:
            # There will always be an Analysis folder
            os.chdir(replicate)
            secondary_dir = os.getcwd()
            print(f"   > Adding {secondary_dir}")

            # Add the header for the first replicate
            with open(charge_file, "r") as current_charge_file:
                if first_charge_file:
                    first_line = current_charge_file.readline()
                    new_charge_file.writelines(first_line)
                    first_charge_file = False
                for index, line in enumerate(current_charge_file):
                    if index == 0:
                        continue
                    elif "nan" in line:
                        print(f"      > Found nan values in {secondary_dir}.")
                    else:
                        new_charge_file.writelines(line)

            os.chdir(primary_dir)

    total_time = round(time.time() - start_time, 3)  # Seconds to run the function
    print(
        f"""
        \t----------------------------ALL RUNS END----------------------------
        \tRESULT: Combined charges across {replicate_count} replicates.
        \tOUTPUT: Generated {charge_file} in the current directory.
        \tTIME: Total execution time: {total_time} seconds.
        \t--------------------------------------------------------------------\n
        """
    )

def get_charge_schemes():
    """
    Calculated a variety of charge schemes with Multiwfn.

    Compute metal-centered Hirshfeld, Voronoi, Mulliken, and ADCH charges.

    Parameters
    ----------

    Returns
    -------

    """
    # Old working directory
    owd = os.getcwd()

    # Populate with names of folders of DFT calculations that were sent to job manager
    list_of_file = ['YIDLOP']
    charge_schemes = {"Hirshfeld": "1", "Voronoi": "2", "Mulliken": "5", "ADCH": "11"}

    # Store extracted multiwfn parameters
    temp_dict = {}

    for f in list_of_file:
        os.chdir(f"{f}/scr")
        subprocess.call("module load multiwfn/noGUI", shell=True)
        command_Z = "module load multiwfn/noGUI"
        command_A = f"/opt/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn {f}.molden"
        final_dest = "/home/manets12/AuNanocage/solvated_guests/charge_files/"

        for key in charge_schemes:
            print(f"Current Key: {key}")
            proc = subprocess.Popen(command_A, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
            calc_command = charge_schemes[key]

            # For atomic charge type corresponding to dict key
            commands = ["7", calc_command, "1", "y", "0", "q"]
            output = proc.communicate("\n".join(commands).encode())
            lines = str(output[0]).split("\\n")
            start = False
            counter = 0

            for num, line in enumerate(lines):
                print(line)
                if "Total valences and free valences" in line:
                    start = True
                    continue
                if start and "Fe" in line:
                    temp_dict[d] = {"Total valence": float(line.split()[-2]),
                                "Free valence": float(line.split()[-1])}
                    counter += 1
            
            # Fename the file as the new, desired name
            new_name = f"{f}_{key}.txt"
            os.rename(f"{f}.chg", new_name)
            shutil.copy(new_name, f"{final_dest}{new_name}")   
        os.chdir(owd)


def calculate_esp():
    """
    Calculate the electrostatic potential (ESP).

    Takes the output from a Multiwfn charge calcualtion and calculates
    the ESP.

    Parameters
    ----------

    Returns
    -------

    """
    df = pd.read_csv("gas_phaseYIDLOP_Voronoi.txt", sep="\s+", names=["Atom", "x", "y", "z", "charge"])
    # Coulombic constant in kg*m**3/(s**4*A**2)
    k = 8.987551*(10**9)

    # Convert each column to list for quicker indexing
    atoms = df["Atom"]
    charges = df["charge"]
    xs = df["x"]
    ys = df["y"]
    zs = df["z"]

    # Pick the index of the atom at which the esp should be calculated
    idx_atom = 0

    # Determine position and charge of the target atom
    xo = xs[idx_atom]
    yo = ys[idx_atom]
    zo = zs[idx_atom]
    chargeo = charges[0]
    total_esp = 0

    # Unit conversion
    A_to_m = 10**(-10)
    KJ_J = 10**-3
    faraday = 23.06   #kcal/(mol*V)
    C_e = 1.6023*(10**-19)
    one_mol = 6.02*(10**23)
    cal_J = 4.184

    for idx in range(0, len(atoms)):
        if idx == idx_atom:
            continue
        else:
            # Calculate esp and convert to units (A to m)
            r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
            total_esp = total_esp + (charges[idx]/r)

    final_esp = k*total_esp*((C_e))*cal_J*faraday   # Note that cal/kcal * kJ/J gives 1
    print(f"{str(final_esp)} kJ/(mol*e)")


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    average_by_residues("all_charges.xls", "template.pdb")
