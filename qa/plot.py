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
from scipy.stats import pearsonr


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


def get_charge_distributions(charge_df, out_file, res_x, res_y, ext, axes_range) -> None:
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

    # Extract the data from the dataframe
    x = charge_df[charge_df.columns[0]]
    y = charge_df[charge_df.columns[1]]

    # Calculate the Pearson's correlation
    corr, _ = pearsonr(x, y)

    # Create the plot
    fig, ax = plt.subplots()
    ax.hist2d(x, y, bins=30, range=axes_range)
    ax.set_aspect("equal")
    
    # Add the Pearson's correlation to the plot
    plt.text(0.55, 0.9, f"Pearson's r: {round(corr,2)}", transform=ax.transAxes, color='white')
    
    plt.xlabel(res_x, fontweight="bold")
    plt.ylabel(res_y, fontweight="bold")
    plt.savefig(f"Analysis/3_coupling/{out_file}", bbox_inches="tight", format=ext, dpi=300)


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
    color = ["#CA4B3F", "#317AB4", "g", "k", "m", "c"]

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
    ext = "svg"
    plt.savefig(f"feature_importance.{ext}", bbox_inches="tight", format=ext)
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

    for index, scheme in enumerate(schemes):
        # Generate a dataframe from the data
        df = pd.read_csv(f"{scheme}.csv")

        # Compute the means and standard deviations
        means = df.mean()
        stds = df.std()

        # create the bar chart with error bars
        plt.bar(
            means.index,
            means.values,
            yerr=stds.values,
            color=colors[index],
            align="center",
            ecolor="black",
            capsize=10,
        )
        plt.xlabel("components", weight="bold")
        plt.ylabel(f"{scheme} ESP kJ/(mol x e)", weight="bold")
        plt.axhline(y=0, color="black", linestyle="--")
        ext = "png"
        plt.savefig(f"{scheme}.{ext}", bbox_inches="tight", format=ext, dpi=300)
        plt.close()


def esp_combined_barchart(schemes, width=5, height=4.5) -> None:
    """
    Plots the ESP for all charge schemes as a barchart with error bars.

    The ESP analysis outputs multicolumn dataframes for components.
    This allows us to compare metal-centered ESP contributions for components.
    This combined version also allows us to compare different charge schemes.

    Paramers
    --------
    schemes : List of strings
        A list of the names of the files that you would like to plot.

    Notes
    -----
    After running qa.manage.collect_esp_components(first, last, step),
    csv files for each charge scheme can be generate in the replicate directory.
    Those files can then be moved to another directory and run this script.

    """
    # Set the width of each bar
    bar_width = 0.2
    
    # Read the CSV files and store them in a list of dataframes
    dfs = [pd.read_csv(scheme) for scheme in schemes]
    
    # Define a list of colors for the bars
    qa.plot.format_plot()
    colors = ["blue", "red", "green", "gray", "cyan", "magenta", "yellow", "orange"]
    
    # Create a list of x-coordinates for each bar
    r_values = [range(len(dfs[0].columns))]
    for i in range(1, len(dfs)):
        r_values.append([x + bar_width for x in r_values[i-1]])
    
    # Plot the bars with error bars
    plt.figure(figsize=(width, height))
    for i, df in enumerate(dfs):
        plt.bar(
            r_values[i],
            df.mean().values,
            yerr=df.std().values,
            width=bar_width,
            color=colors[i % len(colors)],
            align="center",
            ecolor="black",
            capsize=10,
            label=schemes[i].split('.')[0]  # Using filename without extension as label
        )
    
    plt.xlabel("components", weight="bold")
    plt.ylabel("ESP kJ/(mol x e)", weight="bold")
    plt.axhline(y=0, color="black", linestyle="--")
    plt.xticks([r + bar_width * (len(dfs) / 2) for r in r_values[0]], dfs[0].mean().index)
    plt.legend(bbox_to_anchor=(0.6, 0.35), frameon=False)  # raise the legend by 20 pixels
    
    extensions = ["png", "svg"]
    for ext in extensions:
        plt.savefig(f"combined.{ext}", bbox_inches="tight", dpi=300, format=ext)


def plot_rmsd(rmsd_list, labels):
    """
    Plots a bar plot of the RMSD for each analog with its standard error.

    Parameters
    ----------
    rmsd_list: List[List[float]]
        List of lists where each list represents an analog,
        and each list contains RMSDs for each frame.
    labels: List[str]
        List of labels for the x-axis, corresponding to the analog names.

    See Also
    --------
    qa.analyze.compute_rmsd()
    qa.analyze.get_rmsd()

    """
    format_plot()

    # Calculate the mean and standard error for each analog
    rmsd_mean = [np.mean(analog_rmsd) for analog_rmsd in rmsd_list]
    rmsd_std_dev = [np.std(analog_rmsd, ddof=1) for analog_rmsd in rmsd_list]

    # Generate a bar plot for each analog
    x_pos = np.arange(len(rmsd_list))
    plt.bar(x_pos, rmsd_mean, yerr=rmsd_std_dev, align='center', alpha=0.5, capsize=10, color="grey")

    # Set the x-axis labels
    plt.xticks(x_pos, labels)

    # Set the plot labels and title
    plt.ylabel('RMSD (Å)', weight="bold")

    # Save output as a 300 dpi PNG
    ext = "png"
    plt.savefig(f'rmsd_plot.{ext}', bbox_inches="tight", dpi=300)


def time_coupling_plot(charge_df, out_file, res_x, res_y, ext) -> None:
    """
    Compares the charge fluctuations of two residues against time

    Parameters
    ----------
    charge_df : pd.DataFrame
        A dataframe with two columns, each corresponding to a residue.

    See Also
    --------
    qa.analyze.td_coupling()

    """
    # Apply Kulik plotting format
    format_plot()

    # Extract the data from the dataframe and compute the rolling mean
    residue_1 = charge_df[charge_df.columns[0]].rolling(window=5).mean()
    residue_2 = charge_df[charge_df.columns[1]].rolling(window=5).mean()
    
    # Calculate the mean and subtract it from the residue values to center the traces at zero
    residue_1_deviation = residue_1 - residue_1.mean()
    residue_2_deviation = residue_2 - residue_2.mean()
    
    # Conversion factor: 200 single points = 1 picosecond
    frame_to_ps = 1 / 20
    x = np.arange(len(charge_df.index)) * frame_to_ps

    # Create the plot
    plt.plot(x, residue_1_deviation, label=f"{res_x}", color='r')
    plt.plot(x, residue_2_deviation, label=f"{res_y}", color='b')
    plt.xlabel("time (ps)", fontweight="bold")
    plt.ylabel("charge deviation from the mean", fontweight="bold")
    plt.legend()
    plt.savefig(f"Analysis/4_time_coupling/{out_file}", bbox_inches="tight", format=ext, dpi=300)


def esp_dist_plot(esp_choice, xlim=None, ylim=None):
    """
    Creates a scatter plot for a distance and the corresponding ESP.

    Parameters
    ----------
    esp_choice: int
        Column index for the esp_data to be plotted.

    xlim: tuple, optional
        Limits for x-axis in the form (xmin, xmax). If not provided, defaults to None.

    ylim: tuple, optional
        Limits for y-axis in the form (ymin, ymax). If not provided, defaults to None.
    """
    # Load in data
    esp_df = pd.read_csv("Hirshfeld_esp.csv")
    esp_data = esp_df.iloc[:, esp_choice].values
    dist_df = pd.read_csv("centroid_distance.csv")
    dist_data = dist_df.iloc[:, 6].values

    # Create color array
    color_indices = np.repeat(np.arange(len(esp_data) // 400), 400)
    color_indices = np.pad(color_indices, (0, len(esp_data) - len(color_indices)), mode='edge')
    colors = plt.cm.viridis(color_indices / color_indices.max())

    # Create the scatter plot
    format_plot()
    plt.figure(figsize=(5, 5))
    scatter = plt.scatter(dist_data, esp_data, c=colors)
    plt.xlabel("Arg27···Fe Distance (Å)", fontweight="bold")
    plt.ylabel("ESP (kJ/(mol x e))", fontweight="bold")
    
    # Set the x and y limits if specified
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    # Set the xticks to show every integer within xlim
    if xlim is not None:
        plt.xticks(np.arange(np.ceil(xlim[0]), np.floor(xlim[1])+1, 1))
    # If you want even numbers only
        # start = np.ceil(xlim[0]) if np.ceil(xlim[0]) % 2 == 0 else np.ceil(xlim[0]) + 1
        # plt.xticks(np.arange(start, np.floor(xlim[1])+1, 2))

    extensions = ["svg","png"]
    for ext in extensions:
        plt.savefig(f"esp_dist.{ext}", bbox_inches="tight", format=ext, dpi=300)


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    esp_dist_plot()
