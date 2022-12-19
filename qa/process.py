"""Process the raw TeraChem output for future analysis."""

import os
import re
import sys
import glob
import time
import shutil
import lib
from biopandas.pdb import PandasPdb


def get_pdb() -> str:
    """
    Searches all directories recursively for a PDB file.

    If more than one PDB is found it will use the first one.
    If no PDB file was found it will prompt the user for the PDB file path.

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
    pdb_files = sorted(glob.glob("./**/*.pdb", recursive=True))
    for index, pdb in enumerate(pdb_files):
        # Trajectory PDB's should be marked as ensemble or traj
        if "ensemble" in pdb or "traj" in pdb or "top" in pdb:
            continue
        else:
            pdb_file = pdb_files[index]

    # Multiple or no PDB files found scenarios
    if len(pdb_files) == 1:
        print(f"Using {pdb_file} as the template PDB.")
    elif len(pdb_files) > 1:
        print(f"More than one PDB file found -> Using {pdb_file}.")
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

    # Get the xyz from the current directory
    xyz_names = sorted(glob.glob("*.xyz"))

    # Check the results to confirm that there was only one xyz file found
    if len(xyz_names) == 1:
        xyz_name = xyz_names[0]
        print(f"Found the coordinate file {xyz_name}.")
    elif len(xyz_names) > 1:
        xyz_name = xyz_names[0]
        print(f"Found more than one XYZ, using {xyz_name}.")
    else:
        xyz_name = input("No XYZ was found. What is the path to your XYZ? ")

    return xyz_name

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

    # Get the xyz from the current directory
    charge_files = sorted(glob.glob("*.xls"))

    # Check the results to confirm that there was only one .xls file found
    if len(charge_files) == 1:
        charge_file = charge_files[0]
        print(f"Found the charge file {charge_file}.")
    elif len(charge_files) > 1:
        charge_file = charge_files[0]
        print(f"Found more than one .xls, using {charge_file}.")
    else:
        charge_file = input("No .xls was found. What is the path to your .xls? ")

    return charge_file


def combine_runs(
    all_charges: str = "all_charges.xls",
    all_coors: str = "all_coors.xyz",
    atom_count=493,
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
        The number of atoms +2 to account for the header lines.

    Notes
    -----
    Run from the directory that contains the run fragments.
    The variable atom_count should be updated to reflect your system.
    Could be automated in the future by adding to get_PDB.

    """

    # The files where all the charges and coors will be combined
    start_time = time.time()  # Used to report the executation speed
    frame_count = -1  # Report to user with -1 to account for header

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
            for coor in range(0, run_info[index + 1][0] * atom_count):
                all_coors_file.write(coors_file[coor])
        elif (index < len(run_info) - 1) and index != 0:
            for charge in range(1, run_info[index + 1][0] - run[0] + 1):
                all_charges_file.write(charge_file[charge])
                frame_count += 1
            for coor in range(
                0, (run_info[index + 1][0] - run_info[index][0]) * atom_count
            ):
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
    Run all_runs() first for if each replicated was run across multiple runs.
    Generalized to combine any number of replicates.

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
    for (file_name, file_location) in zip(new_file_names, file_locations):
        # Open a new file where we will right the concatonated output
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


def xyz2pdb_traj() -> None:
    """
    Converts an xyz trajectory file into a pdb trajectory file.

    Note
    ----
    Assumes that the only TER flag is at the end.

    """

    # Search for the XYZ and PDB files names
    start_time = time.time()  # Used to report the executation speed
    pdb_name = get_pdb()
    xyz_name = get_xyz()
    # Remove the extension to get the protein name to use as the PDB header
    protein_name = pdb_name.split("/")[-1][:-4]
    new_pdb_name = f"{protein_name}_traj.pdb"

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
                print(f"Script died at {line_count} -> '{line}'")
                quit()
            pdb_line = pdb_file[atom - 2]  # PDB is two behind the xyz
            new_file.write(
                f"{pdb_line[0:30]}{x[0:6]}  {y[0:6]}  {z[0:6]}  {pdb_line[54:80]}\n"
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
    Converts an xyz trajectory file into a pdb trajectory file.

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


def remove_incomplete_xyz() -> None:
    """
    For removing incomplete frames during troublshooting.

    It seems that sometimes an frame was not written out incompletely.
    In these cases, the best course is to remove them.

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


def check_valid_resname(res) -> tuple[str, int]:
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
    aa_nume : int
        The requested amino acid's position in the sequence.
    
    """
    # Get amino acid identifier information from library module
    aa_identifiers = lib.get_aa_identifiers()
    # Check if one or three letter code provided by counting letters
    letter_count = sum(map(str.isalpha, res))

    if letter_count == 1:
        aa_name = res[0].upper().strip() # clean the input
        aa_num = int(res[1:])
        # Check if the provided one-letter code matches a known amino acid
        if not any([True for k,v in aa_identifiers.items() if v[0] == aa_name]):
            raise ValueError("The provided one-letter code is unknown.")
    
    if letter_count == 3:
        aa_name = res[:3].upper().strip() # clean the input
        aa_num = int(res[3:])
        # Check if the provided three-letter code matches a known amino acid
        if not any([True for k,v in aa_identifiers.items() if v[1] == aa_name]):
            raise ValueError("The provided three-letter code is unknown.")

    print(f"Reqesting amino acid {aa_name} at index {aa_num}.")

    return aa_name, aa_num


def get_res_atom_indices(res, scheme="all") -> list[int]:
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
    residue_df = ppdb[(ppdb["residue_name"] == aa_name) & (ppdb["residue_number"] == aa_num)]
    atom_index_list = residue_df.index.tolist()
    print(residue_df)
    print(atom_index_list)
    # Use if you only want the backbone atoms summed
    if scheme == "backbone":
        print("Retrieving only backbone indices.")
        bb_atoms = ["N", "H", "C", "O"]
        backbone_df = residue_df[residue_df["atom_name"].isin(bb_atoms)]
        atom_index_list = backbone_df.index.tolist()
    if scheme != "all" and scheme != "backbone":
        raise ValueError("Scheme not recognized. Select all or backbone.")
    # Alert the user if the list comes out empty
    if len(atom_index_list) == 0:
        raise ValueError("ERROR: No atom indices were found. Verify that it exists.")

    return atom_index_list


if __name__ == "__main__":
    # Run when executed as a script
    get_res_atom_indices("Lys8")
