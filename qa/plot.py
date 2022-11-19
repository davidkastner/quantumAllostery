"""Prebuilt vizualization functions."""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np


def format_plot():
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


def heatmap(csv):
    """
    Generates formatted heat maps.

    Uses the Kulik Lab figure formatting standards.

    """

    # General styling variables
    light_gray = "#8e8e8e"
    dark_gray = "#7e7e7e"
    residues = [
        "ACE1",
        "ASP2",
        "GLU3",
        "GLN4",
        "GLN5",
        "LEU6",
        "HIS7",
        "SER8",
        "GLN9",
        "LYS10",
        "ARG11",
        "LYS12",
        "ILE13",
        "THR14",
        "LEU15",
        "NHE16",
        "ACE17",
        "ASP18",
        "GLU19",
        "GLN20",
        "GLN21",
        "LEU22",
        "SER23",
        "SER24",
        "GLN25",
        "LYS26",
        "ARG27",
        "NHE28",
        "HEME",
        "FE",
    ]

    # Generate dataset
    matrix = np.genfromtxt(csv, delimiter=",")
    np.fill_diagonal(matrix, 0)
    df = pd.DataFrame(matrix)
    df.columns = residues
    df.index = residues

    # Apply base Kulik plot parameters
    format_plot()

    # Generate plot
    ax = sns.heatmap(
        df,
        cmap="RdBu",
        vmin=-0.5,
        vmax=0.5,
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
    ax.hlines([0, 30], colors="k", *ax.get_xlim(), linewidth=3.5)
    ax.vlines([0, 30], colors="k", *ax.get_ylim(), linewidth=3.5)

    # plt.savefig("final.png", bbox_inches="tight", format="png", dpi=300)
    plt.savefig("final.svg", bbox_inches="tight", format="svg")


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    heatmap("chargematbb.csv")
