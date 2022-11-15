"""Provide the primary functions."""

import process
import predict
import plot


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


def process():
    """
    Placeholder function to show example docstring (NumPy format).

    Replace this function and doc string for your own project.

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from.

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution.
    """

    return


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
