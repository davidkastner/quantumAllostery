"""Command-line interface entry point."""

# Print first to welcome the user while it waits to load the modules
print("\n.-----------------------------------.")
print("| WELCOME TO QUANTUM ALLOSTERY (QA) |")
print(".-----------------------------------.")
print("Default programmed actions for the quantumAllostery package.")
print("Feel free to add your own action sequences.")
print("GitHub: https://github.com/davidkastner/quantumAllostery")
print("Documenation: https://quantumallostery.readthedocs.io\n")
print("Loading...")

# !pip show biopandas

import sys
import qa.analyze
import qa.library
import qa.predict
import qa.process
import qa.plot
import qa.logger

def cli():
    """
    The overall command-line interface (CLI) entry point.
    The CLI interacts with the rest of the package.

    It will prompt the user with multiple levels of options.
    This is advantagous as quantumAllostery is moderate in scope,
    and because it introduces the user to the available functionality.
    Improves long-term maintainability.

    Notes
    -----
    The functions will be expecting the following file structure:
    .
    ├── variant
    │   ├── replicate
    │   │   ├── start
    │   │   ├── restart
    │   │   └── ...
    │   ├── replicate
    │   │   ├── start
    │   │   ├── restart
    │   │   └── ...
    │   └── ...
    ├── variant
    │   ├── replicate
    │   │   ├── start
    │   │   ├── restart
    │   │   └── ...
    │   ├── replicate
    │   │   ├── start
    │   │   ├── restart
    │   │   └── ...
    │   └── ...
    └── ...

    """

    # First level of actions
    print("""
        1) Process files
        2) Analyze results
        3) Predict using ML model
        4) Plot results
        5) Help
        6) Quit
        """)
    request = input("Select an option: ")

    # Second level of actions
    if request == "1": #Process files
        print("""
        1) Combine restarted trajectories
        2) Combine restarted trajectories for multiple replicates
        3) Combine replicates into a single trajectory
        4) Convert an XYZ to PDB trajectory
        5) Remove incomplete xyz frames
        6) Help
        7) Quit
        """)
        request = input("Select an option: ")
        
        if request == "1": # Combine restarted trajectories
            atom_count = qa.process.get_atom_count()
            qa.process.combine_restarts(atom_count)
        elif request == "2": # Combine restarted trajectories for multiple replicates
            atom_count = qa.process.get_atom_count()
            qa.process.run_all_replicates(lambda: qa.process.combine_restarts(atom_count))
        elif request == "3": # Combine replicates into a single trajectory
            qa.process.combine_replicates()
        elif request == "4": # Convert an XYZ to PDB trajectory
            qa.process.xyz2pdb_traj()
        elif request == "5": # Remove incomplete xyz frames
            qa.process.remove_incomplete_xyz()
        elif request == "6": # Help
            qa.logger.help()
        elif request == "7": # Quit
            qa.logger.clean_exit()
        else:
            qa.logger.nonoption_exit()

    elif request == "2": # Analyze results
        print("""
        1) Quit
        """)
        request = input("Select an option: ")

    elif request == "3": # Predict using ML model
        print("""
        1) Quit
        """)
        request = input("Select an option: ")

    elif request == "4": # Plot results
        print("""
        1) Plot charge-coupling heatmap for two residues
        2) Quit
        """)
        request = input("Select an option: ")

        if request == "1": # Combine restarted trajectories
            res_x = input("\n> What is the first residue (Asp1)? ")
            res_y = input("> What is the second residue (Gly2)? ")
            charge_df = qa.analyze.get_joint_qres(res_x, res_y)
        elif request == "2": # Quit
            qa.logger.clean_exit()
        else:
            qa.logger.nonoption_exit()

    elif request == "1": # Plot charge-coupling heatmap for two residues
        qa.logger.help()
    elif request == "2": # Quit
        qa.logger.clean_exit()
    else: # Non-option
        qa.logger.nonoption_exit()


def level_one():
    """
    Lists all the potential modules to choose from.
    """
        # First level of actions
    print("""
        1) Process files
        2) Analyze results
        3) Predict using ML model
        4) Plot results
        5) Help
        6) Quit
        """)
    request = input("Select an option: ")

if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
