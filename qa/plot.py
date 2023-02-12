"""Prebuilt vizualization functions."""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from typing import List
from matplotlib.ticker import MultipleLocator
import numpy as np
import qa.library


def format_plot() -> None:
    """
    General plotting parameters for the Kulik Lab.

    """
    font = {"family": "sans-serif", "weight": "bold", "size": 10}
    plt.rc("font", **font)
    plt.rcParams['xtick.major.pad'] = 5
    plt.rcParams['ytick.major.pad'] = 5
    plt.rcParams["axes.linewidth"] = 2
    plt.rcParams["xtick.major.size"] = 7
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 7
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["svg.fonttype"] = "none"


# def heatmap(data, protein, delete, out_file, cmap) -> None:
def heatmap(data: str, residues: List[str], delete: List[List[int]], out_file: str="heatmap.svg", cmap="RdBu", v=[-.4,.4]) -> None:
    """
    Generates formatted heat maps.

    Uses the Kulik Lab figure formatting standards.

    Parameters
    ----------
    data: str
        The name of the data file, which can be a data or a dat file.
    delete: List[List[int]]
        A list of the amino acids you would like removed indexed at zero.
    out_file: str
        The name you would like the image saved as.
    cmp: str
        Your color map setting.

    Examples
    --------
    heatmap(data="cacovar.dat", protein="mc6sa", delete=[0,15,16,27,28,29], out_file="matrix_geom.svg")

    """

    # General styling variables
    light_gray = "#8e8e8e"
    dark_gray = "#7e7e7e"

    # Delete extra residues
    residues = [item for index, item in enumerate(residues) if index not in delete[0]]

    # Identify matrix format and read in
    contents = open(data, "r").readlines()
    contents_joined = " ".join(contents)  # Create a single string for parsing
    if "," in contents_joined: # If a csv file
        matrix = np.genfromtxt(data, delimiter=",")
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
        df = df.drop(df.columns[delete[1]], axis=0)
        df = df.drop(df.columns[delete[1]], axis=1)
    df.columns = residues
    df.index = residues

    # Apply base Kulik plot parameters
    format_plot()

    # Generate plot
    ax = sns.heatmap(
        df,
        cmap=cmap,
        vmin=v[0],
        vmax=v[1],
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


def get_parity_plot(x: List[int], y: List[int]) -> None:
    """
    General set up to create an attractive parity plot.

    Parameters
    ----------
    x: list[int]
        List of x data points.
    y: list[int]
        List of y data points.

    """

    # User defined values
    format_plot()
    xlabel = "UB3LYP/def2-TZVP"
    ylabel = "UB3LYP/LACVP*"

    # This is just example data, replace it with your own

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_aspect("equal", adjustable="box")  # Parity plots should be square
    # Create the identity line
    plt.plot([-2, 40], [-2, 40], color="r", linewidth=2, zorder=0)
    # Create scatter plot of data
    plt.scatter(x, y, marker="o", color="k", zorder=1)
    plt.xlabel(xlabel, fontsize=10, weight="bold")
    plt.ylabel(ylabel, fontsize=10, weight="bold")
    plt.tick_params(which="both", bottom=True, top=True, left=True, right=True)
    # Force the axes to be the same
    plt.xticks(np.arange(0, 36, 5))
    plt.yticks(np.arange(0, 36, 5))
    ax.set_xlim(0, 37)
    ax.set_ylim(0, 37)
    # Create the plot
    plt.savefig(
        "parity.png", format="png", dpi=600, bbox_inches="tight", transparent=True
    )


def get_charge_distributions(charge_df, out_file, res_x, res_y, ext):
    """
    Creates a charge distribution plot with one residue on each axis.

    Parameters
    ----------
    charge_pd : pd.DataFrame
        A dataframe with two columns, each corresponding to a residue.

    See Also
    --------
    qa.analyze.get_joint_qres()

    """

    # Apply Kulik plotting format
    format_plot()

    # Define the bins automatically (for testing)
    # x_min = charge_df[charge_df.columns[0]].min()
    # x_max = charge_df[charge_df.columns[0]].max()
    # y_min = charge_df[charge_df.columns[1]].min()
    # y_max = charge_df[charge_df.columns[1]].max()
    # x_bins = np.linspace(x_min, x_max, 50)
    # y_bins = np.linspace(y_min, y_max, 50)

    # Extract the data from the dataframe
    x = charge_df[charge_df.columns[0]]
    y = charge_df[charge_df.columns[1]]

    # Create the plot
    fig, ax = plt.subplots()
    # ax.hist2d(x, y, bins=(x_bins, y_bins))
    ax.hist2d(x, y, bins=30, range=[[0.4, 1.2], [-1.2, -0.4]])
    ax.set_aspect('equal')
    plt.xlabel(res_x, fontweight='bold')
    plt.ylabel(res_y, fontweight='bold')
    plt.savefig(out_file, bbox_inches="tight", format=ext, dpi=300)
