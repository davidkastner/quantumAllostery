"""Command-line interface entry point."""

import sys
import qa.process
import qa.predict
import qa.plot


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

    # User welcome header
    print("\n.-----------------------------------.")
    print("| WELCOME TO QUANTUM ALLOSTERY (QA) |")
    print(".-----------------------------------.\n")
    print("Default programmed actions for the quantumAllostery package.")
    print("Feel free to add your own action sequences.")
    print("\t• GitHub: https://github.com/davidkastner/quantumAllostery")
    print("\t• Documenation: https://quantumallostery.readthedocs.io/en/latest/\n")

    # First level of actions
    lvl_1 = """
        1) Process raw files
        2) Predict charge-coupling interactions
        3) Visualize results
        4) Quit\n
        """
    print(f"{lvl_1}")
    request = input("Select an option: ")

    # Second level of actions
    if request == "1":
        lvl_11 = """
        1) Combine restarted trajectories
        2) Combine replicates into a single trajectory
        3) Convert an XYZ to PDB trajectory
        4) Convert an XYZ to PDB ensemble
        5) Remove incomplete xyz frames
        6) Quit\n
        """
        print(f"{lvl_11}")
        request = input("Select an option: ")
        
        if request == "1":
            qa.process.combine_restarts()
        elif request == "2":
            qa.process.combine_replicates()
        elif request == "3":
            qa.process.xyz2pdb_traj()
        elif request == "4":
            qa.process.xyz2pdb_ensemble()
        elif request == "5":
            qa.process.remove_incomplete_xyz()
        elif request == "6":
            sys.exit("Exited successfully.\n Thanks for using Quantum Allostery.\n")
        else:
            print("Sorry that option is not available.")
            print("Please select one of the available options.\n")

    elif request == "2":
        lvl_12 = """
        1) Quit\n
        """
        print(f"{lvl_12}")
        request = input("Select an option: ")

    elif request == "3":
        lvl_13 = """
        1) Process raw files
        2) Predict charge-coupling interactions
        3) Visualize results
        4) Quit\n
        """
        print(f"{lvl_13}")
        request = input("Select an option: ")

    elif request == "3":
        sys.exit("Exited successfully.\n Thanks for using Quantum Allostery.\n")

    else:
        print("Sorry that option is not available.")
        print("Please select one of the available options.\n")


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
