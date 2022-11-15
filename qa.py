"""Command-line interface entry point."""

import qa.process as process
import qa.predict as predict
import qa.plot as plot


def cli():
    """
    The overall command-line interface (CLI) that interacts with the rest of the package.

    """
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
