"""Prebuilt vizualization functions."""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np


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
    plt.rcParams['svg.fonttype'] = 'none'


# def heatmap(csv: str, protein: str, delete: list[int]=[], out_file: str="heatmap.svg", cmap="RdBu") -> None:
def heatmap(csv, protein, delete, out_file, cmap) -> None:
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
    contents_joined = " ".join(contents)  # Create a single string for parsing
    if "," in contents_joined:  # If CSV
        matrix = np.genfromtxt(csv, delimiter=",")
    elif " " in contents_joined:  # If a dat file
        matrix = []
        for line in contents:
            matrix.append(line.split())
        matrix = [[float(j) for j in i] for i in matrix]  # Strings to float
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

    ext = out_file.split(".")[-1]  # Determine the output file type
    plt.savefig(out_file, bbox_inches="tight", format=ext, dpi=300)
    
def get_parity_plot():
    """
    General set up to create an attractive parity plot.
    """
    # User defined values
    format_plot()
    xlabel = 'UB3LYP/def2-TZVP'
    ylabel = 'UB3LYP/LACVP*'
    x = [0.00,7.39,20.62,0.00,3.16,25.82,0.00,14.37,17.15,0.00,5.08,17.78,0.00,3.88,32.40,0.00,6.85,35.26,0.00,6.14,27.07,0.00,-1.06,33.19]
    y = [0.00,7.63,21.88,0.00,6.79,25.05,0.00,14.50,17.83,0.00,6.73,18.80,0.00,2.77,30.52,0.00,7.06,31.66,0.00,5.84,25.29,0.00,0.15,30.92]

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_aspect('equal', adjustable='box') # Parity plots should be square
    # Create the identity line
    plt.plot([-2,40],[-2,40], color="r", linewidth=2, zorder=0)
    # Create scatter plot of data
    plt.scatter(x, y, marker="o", color="k", zorder=1)
    plt.xlabel(xlabel, fontsize=10, weight="bold")
    plt.ylabel(ylabel, fontsize=10, weight="bold")
    plt.tick_params(which="both", bottom=True, top=True, left=True, right=True)
    # Force the axes to be the same
    plt.xticks(np.arange(0, 36, 5))
    plt.yticks(np.arange(0, 36, 5))
    ax.set_xlim(0,37)
    ax.set_ylim(0,37)
    # Create the plot
    plt.savefig("parity.png", format="png", dpi=600, bbox_inches='tight', transparent=True)


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    #  delete = [0,15,16,27,28,29]
    heatmap(csv="cacovar.dat", protein="mc6sa", out_file="matrix_geom.svg")
