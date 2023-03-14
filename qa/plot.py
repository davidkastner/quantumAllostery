"""Prebuilt vizualization functions."""

import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from typing import List
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import numpy as np
import qa.process


def format_plot() -> None:
    """
    General plotting parameters for the Kulik Lab.

    """
    font = {"family": "sans-serif", "weight": "bold", "size": 10}
    plt.rc("font", **font)
    plt.rcParams["xtick.major.pad"] = 5
    plt.rcParams["ytick.major.pad"] = 5
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
def heatmap(
    data: str,
    residues: List[str],
    delete: List[List[int]],
    out_file: str = "heatmap.svg",
    cmap="RdBu",
    v=[-0.4, 0.4],
) -> None:
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
    if "," in contents_joined:  # If a csv file
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
    plt.close()


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


def get_charge_distributions(charge_df, out_file, res_x, res_y, ext) -> None:
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

    # # Define the bins automatically (for testing)
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
    ax.hist2d(x, y, bins=30, range=[[-0.2, 0.44], [-1.25, -0.53]])
    ax.set_aspect("equal")
    plt.xlabel(res_x, fontweight="bold")
    plt.ylabel(res_y, fontweight="bold")
    plt.savefig(
        f"Analysis/3_coupling/{out_file}", bbox_inches="tight", format=ext, dpi=300
    )
    plt.close()


def plot_feature_importance(
    models: List[str], template: List[str], mutations: List[str], by_atom=False
) -> None:
    """
    Creates a plot with the features on the x-axis and their importance on the y

    Parameters
    ----------
    models: List[str]
        A list of the different models that were trained and tested.
    template: str
        The name of the template pdb for the protein of interest.

    See Also
    --------
    qa.predict.run_ml()

    """
    # Apply Kulik plotting format
    format_plot()

    # Get the amino acid names of our features
    residues_indentifier = qa.process.get_residue_identifiers(template, by_atom=by_atom)
    residues_indentifier = [
        x for i, x in enumerate(residues_indentifier) if i not in mutations
    ]

    # Get the feature importance data which has been stored by Demystifying
    root = os.getcwd()
    feature_sets = []
    for model in models:
        feature_importance = list(np.load(f"{root}/{model}/feature_importance.npy"))
        feature_sets.append(feature_importance)

    # Create the plot
    print(f"   > Creating a plot of feature importance for all models.")
    x_axis = residues_indentifier
    color = ["b", "r", "g", "k", "m", "c"]

    # Requesting a by atom plot?
    if by_atom:
        x_axis = [x for x in range(len(feature_sets[0]))]

    # feature_set is a list of list of lists [[[1,2,1], [3,2,2]]]
    averages = [
        [np.mean(sublist) for sublist in list_of_lists]
        for list_of_lists in feature_sets
    ]
    errors = [
        [np.std(sublist) for sublist in list_of_lists] for list_of_lists in feature_sets
    ]

    fig, ax = plt.subplots()

    # Show the average with the error as a transparent trace
    for index, average in enumerate(averages):
        plt.plot(x_axis, average, color=color[index], linewidth=2)
        upper = np.array(average) + np.array(errors[index])
        lower = np.array(average) - np.array(errors[index])
        plt.fill_between(x_axis, lower, upper, color=color[index], alpha=0.1)

    # Add final touches
    xlabel = "residues"
    ylabel = "feature importance score"
    plt.xlabel(xlabel, weight="bold")
    plt.ylabel(ylabel, weight="bold")
    plt.tick_params(axis="x", which="major", rotation=90)
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax.tick_params(axis="y", which="minor", length=4, width=2)
    plt.savefig("feature_importance.png", bbox_inches="tight", format="png", dpi=300)
    plt.close()

def esp_separate_barchart() -> None:
    """
    Plots the ESP for one single charge scheme as a barchart with error bars.

    The ESP analysis outputs multicolumn dataframes for components.
    This allows us to compare metal-centered ESP contributions for components.

    """
    # Apply Kulik plotting format
    qa.plot.format_plot()

    schemes = ["ADCH", "Hirshfeld", "Mulliken", "Voronoi"]
    colors = ["blue", "red", "green", "gray"]

    for index,scheme in enumerate(schemes):
        # Generate a dataframe from the data
        df = pd.read_csv(f"{scheme}.csv")

        # Compute the means and standard deviations
        means = df.mean()
        stds = df.std()
        
        # create the bar chart with error bars
        plt.bar(means.index, means.values, yerr=stds.values, color=colors[index], align='center', ecolor='black', capsize=10)
        plt.xlabel('Components', weight="bold")
        plt.ylabel(f'{scheme} ESP kJ/(mol x e)', weight="bold")
        plt.axhline(y=0, color='black', linestyle='--')
        ext = "png"
        plt.savefig(f"{scheme}.{ext}", bbox_inches="tight", format=ext, dpi=300)
        plt.close()


def esp_combined_barchart() -> None:
    """
    Plots the ESP for all charge schemes as a barchart with error bars.

    The ESP analysis outputs multicolumn dataframes for components.
    This allows us to compare metal-centered ESP contributions for components.
    This combined version also allows us to compare different charge schemes.
    
    """
    # Apply Kulik plotting format
    qa.plot.format_plot()

    # Set the width of each bar
    bar_width = 0.2

    schemes = ["ADCH", "Hirshfeld", "Mulliken", "Voronoi"]
    df1 = pd.read_csv(f'{schemes[0]}.csv')
    df2 = pd.read_csv(f'{schemes[1]}.csv')
    df3 = pd.read_csv(f'{schemes[2]}.csv')
    df4 = pd.read_csv(f'{schemes[3]}.csv')

    
    # Create a list of x-coordinates for each bar
    r1 = range(len(df1.columns))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
    r4 = [x + bar_width for x in r3]

    # create the bar chart with error bars
    plt.bar(r1, df1.mean().values, yerr=df1.std().values, width=bar_width, color="blue", align='center', ecolor='black', capsize=10, label=schemes[0])
    plt.bar(r2, df2.mean().values, yerr=df2.std().values, width=bar_width, color="red", align='center', ecolor='black', capsize=10, label=schemes[1])
    plt.bar(r3, df3.mean().values, yerr=df3.std().values, width=bar_width, color="green", align='center', ecolor='black', capsize=10, label=schemes[2])
    plt.bar(r4, df4.mean().values, yerr=df4.std().values, width=bar_width, color="gray", align='center', ecolor='black', capsize=10, label=schemes[3])
    
    plt.xlabel('Components', weight="bold")
    plt.ylabel('ESP kJ/(mol x e)', weight="bold")
    plt.axhline(y=0, color='black', linestyle='--')
    
    plt.xticks([r + bar_width*1.5 for r in r1], df1.mean().index)
    plt.legend()

    ext = "png"
    plt.savefig(f"combined.{ext}", bbox_inches="tight", format=ext, dpi=300)
    plt.close()

if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    esp_nomulliken_barchart()