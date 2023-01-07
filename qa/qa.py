"""Command-line interface entry point."""

import process
import predict
import plot


def cli():
    """
    The overall command-line interface (CLI) that interacts with the rest of the package.

    Notes
    -----
    Expecting the following file structure:
    .
    ├── variant 1
    │   ├── replicate 1
    │   │   ├── segment 1
    │   │   ├── segment 2
    │   │   └── ...
    │   ├── replicate 2
    │   │   ├── segment 1
    │   │   ├── segment 2
    │   │   └── ...
    │   └── ...
    ├── variant 2
    │   ├── replicate 1
    │   │   ├── segment 1
    │   │   ├── segment 2
    │   │   └── ...
    │   ├── replicate 2
    │   │   ├── segment 1
    │   │   ├── segment 2
    │   │   └── ...
    │   └── ...
    └── ...

    """

    # User welcome header
    print("\n.----------------------------------.")
    print("| WELCOME TO QUANTUM ALLOSTERY (QA)|")
    print(".----------------------------------.\n")
    print("\tGitHub: https://github.com/davidkastner/quantumAllostery")
    print("\tDocumenation: https://quantumallostery.readthedocs.io/en/latest/\n")

    # Ask the user what task to perform with quantumAllostery
    prompt = "What would you like to do (e.g., 1, 2, or 3)?\n"
    option_1 = "\t1) Process raw files\n"
    option_2 = "\t2) Predict charge-coupling interactions\n"
    option_3 = "\t3) Visualize results\n"
    request = input(f"{prompt}{option_1}{option_2}{option_3}")

    # Execute requested task
    if request == "1":
        process()
    elif request == "2":
        predict()
    elif request == "3":
        plot()
    else:
        print("Select one of the available options.")


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
