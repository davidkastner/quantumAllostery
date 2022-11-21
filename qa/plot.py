"""Prebuilt vizualization functions."""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
# import lib


def format_plot() -> None:
    """
    General plotting parameters for the Kulik Lab.

    """
    font = {"family": "sans-serif", "weight": "bold", "size": 10}
    plt.rc("font", **font)
    plt.rcParams["axes.linewidth"] = 2
    plt.rcParams["xtick.major.size"] = 7
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 7
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["xtick.direction"] = "out"
    plt.rcParams["ytick.direction"] = "in"


# def heatmap(csv: str, protein: str, delete: list[int]=[], out_file: str="heatmap.svg", cmap="RdBu") -> None:
def heatmap(csv, protein, delete, out_file, cmap):
    """
    Generates formatted heat maps.

    Uses the Kulik Lab figure formatting standards.

    """
    
    # General styling variables
    light_gray = "#8e8e8e"
    dark_gray = "#7e7e7e"
    remove: list(int) = []
    residues = lib.sequence(protein)

    # Identify matrix format and read in
    contents = open(csv, "r").readlines()
    contents_joined = " ".join(contents) # Create a single string for parsing
    if "," in contents_joined:  # If CSV
        matrix = np.genfromtxt(csv, delimiter=",")
    elif " " in contents_joined:  # If a dat file
        matrix = []
        for line in contents:
            matrix.append(line.split())
        matrix = [[float(j) for j in i] for i in matrix] # Strings to float
        matrix = np.array(matrix)

    np.fill_diagonal(matrix, 0)  # Set the diagonal to zero as they are trivial

    df = pd.DataFrame(matrix)
    # Remove specific rows and columns from non-residues
    if len(delete) > 0:
        df = df.drop(df.columns[delete], axis=0)
        df = df.drop(df.columns[delete], axis=1)
    df.columns = residues
    df.index = residues

    # Apply base Kulik plot parameters
    format_plot()

    # Generate plot
    ax = sns.heatmap(
        df,
        cmap=cmap,
        vmin=-1,
        vmax=1,
        xticklabels=True,
        yticklabels=True,
        linewidth=0.03,
        linecolor=light_gray,
    )
    ax.get_yaxis().set_tick_params(direction="out")
    plt.ylabel("residues", fontsize=10, weight="bold")
    plt.xlabel("residues", fontsize=10, weight="bold")
    plt.tick_params(axis="y", which="major", labelsize=8, rotation=0, length=0)
    plt.tick_params(axis="x", which="major", labelsize=8, rotation=90, length=0)
    ax.xaxis.tick_top()  # x axis on top
    ax.xaxis.set_label_position("top")
    # Add lines every five residues
    ax.hlines([5, 10, 15, 20, 25], colors=dark_gray, *ax.get_xlim(), linewidth=1.5)
    ax.vlines([5, 10, 15, 20, 25], colors=dark_gray, *ax.get_ylim(), linewidth=1.5)
    # Add broders
    ax.hlines([0, len(residues)], colors="k", *ax.get_xlim(), linewidth=3.5)
    ax.vlines([0, len(residues)], colors="k", *ax.get_ylim(), linewidth=3.5)

    ext = out_file.split(".")[-1] # Determine the output file type
    plt.savefig(out_file, bbox_inches="tight", format=ext, dpi=300)


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    #  delete = [0,15,16,27,28,29]
    heatmap(csv="cacovar.dat", protein="mc6sa", out_file="matrix_geom.svg")
