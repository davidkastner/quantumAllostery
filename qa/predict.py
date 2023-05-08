"""Layer-wise relevance propogate MD predictions."""

import os
import sys
import pandas as pd
import numpy as np
from typing import List, Tuple
from demystifying import feature_extraction as fe, visualization
from demystifying import relevance_propagation as relprop
from sklearn.preprocessing import MinMaxScaler
from sklearn.utils import shuffle
import qa.process


def shuffle_data(charges_mat, labels_mat):
    """
    Enables the option to shuffle the time series data if desired.
    """
    # Shuffle the data in groups of 100
    n_samples = charges_mat.shape[0]
    n_samples = int(n_samples / 100) * 100
    inds = np.arange(n_samples)
    inds = inds.reshape((int(n_samples / 100), 100))
    perm_inds = np.random.permutation(inds)
    perm_inds = np.ravel(perm_inds)
    samples = charges_mat[perm_inds]  # Apply the shuffling to the matrix
    labels = labels_mat[perm_inds]  # Apply the same shuffling to the labels
    print("   > Shuffling data in blocks of 100.")

    return samples, labels


def create_combined_csv(
    charge_files: List[str], templates: List[str], mutations: List[int]
) -> pd.DataFrame:
    """
    Generate a pd.DataFrame of all features.

    Returns
    -------
    charges_df: pd.DataFrame
        The original charge data as a pandas dataframe.
    lablels_df: pd.DataFrame
        One-hot-encoded labels for each frame.

    """
    # Convert the input data files to pd.DataFrames
    dataframes = []
    labels_df = pd.DataFrame()
    for charge_file, template in zip(charge_files, templates):
        print(f"   > Converting atoms to residues for {charge_file}.")
        # Average the charges by residue
        # We does this to minimize the inaccuracies of mulliken charges
        avg_by_residues = qa.process.average_charge_residues(charge_file, template)

        # Add a column for the one-hot-encoded labels for each frame
        print(f"   > Creating labels for {charge_file}.")
        label = [1 for x in range(len(avg_by_residues))]
        avg_by_residues[f"{charge_file}"] = label

        # Store the labeled, averaged frames in the list
        dataframes.append(avg_by_residues)

        # Remove the last column (label) before saving
        avg_by_residues_no_label = avg_by_residues.iloc[:, :-1]

        # Save the individual dataframe without the label column to a CSV file
        output_file = charge_file.split(".")[0] + "_charges.csv"
        avg_by_residues_no_label.to_csv(output_file, index=False)
        print(f"   > Saved {output_file}.")

    # Drop the residue columns that were mutated
    # We can't compare these residues' charges as their atom counts differ
    clean_dataframes = []
    for index, df in enumerate(dataframes):
        df = df.drop(df.columns[mutations], axis=1)
        clean_dataframes.append(df)

    # Combine both dataframes into a single dataframe
    print(f"   > Combining {charge_files}.")
    combined_df = pd.concat(clean_dataframes, ignore_index=True, sort=False)
    combined_df.fillna(0, inplace=True)

    # Break off the last n columns (labels) and save them as their own df
    charges_df = combined_df.iloc[:, : -len(charge_files)]
    labels_df = combined_df.iloc[:, -len(charge_files) :]

    return charges_df, labels_df


def data_processing(df, labels, n_frames=1):
    """
    Scales the data for the ML workflows.

    Parameters
    ----------
    df: pd.DataFrame
        The original data as a pandas dataframe
    labels_df: pd.DataFrame
        One-hot-encoded labels for each frame.
    n_frames: int
        Step for filtering the data (e.g. 1 = every frame, 2 = everyother frame)

    Returns
    -------
    df_norm: pd.DataFrame
        The data scaled by column.

    """
    # Scale each column such that all values are between 0 and 1
    scaler = MinMaxScaler()
    df_norm = pd.DataFrame(
        scaler.fit_transform(df.values), columns=df.columns, index=df.index
    )

    # Convert to numpy matrices for compatibility with Demystify
    data_norm = df_norm.to_numpy()
    labels = labels.to_numpy().astype(int)

    # Filter the data as every n frames
    filtered_data = data_norm[::n_frames]
    filtered_labels = labels[::n_frames]

    return filtered_data, filtered_labels


def run_ml(data_norm, labels, models=["RF", "MLP"], recompute=False):
    """
    ML analysis workflow.

    Parameters
    ----------
    df_norm: numpy matrix
        The data scaled by column.
    labels_df: pd.DataFrame
        One-hot-encoded labels for each frame.

    """
    # Check if the models and feature importances have already been created
    output_exists = [os.path.isdir(model) for model in models]
    if not recompute and not False in output_exists:
        print(f"> Results exist for the following models: {models}")
        return
    else:
        # Output has not been generated yet
        print("   > Creating csv's to check.")
        np.savetxt("data.csv", data_norm, delimiter=",")
        np.savetxt("labels.csv", labels, fmt="%i", delimiter=",")

        # Set the arguments for the ML workflows
        kwargs = {
            "samples": data_norm,
            "labels": labels,
            "filter_by_distance_cutoff": False,
            "lower_bound_distance_cutoff": 1.0,
            "upper_bound_distance_cutoff": 1.0,
            "use_inverse_distances": False,
            "n_splits": 3,
            "n_iterations": 1,
            "scaling": True,
        }

        # Running various ML workflows
        models = models
        feature_extractors = [
            fe.RandomForestFeatureExtractor(
                one_vs_rest=True, classifier_kwargs={"n_estimators": 100}, **kwargs
            ),
            fe.MlpFeatureExtractor(
                classifier_kwargs={
                    "hidden_layer_sizes": (120,),
                    "solver": "adam",
                    "max_iter": 1000000,
                },
                one_vs_rest=True,
                activation=relprop.relu,
                **kwargs,
            ),
        ]

        # Process the results
        postprocessors = []
        working_dir = "."
        for extractor, model in zip(feature_extractors, models):
            print(f"> Running {model} model.")
            extractor.extract_features()
            # Post-process data (rescale and filter feature importances)
            postprocessor = extractor.postprocessing(
                working_dir=working_dir,
                rescale_results=True,
                filter_results=False,
                feature_to_resids=None,
            )
            postprocessor.average()
            postprocessor.evaluate_performance()
            postprocessor.persist()
            postprocessors.append(postprocessor)


# Execute the script
if __name__ == "__main__":
    # Which amino acids to remove
    mutations = [2,19,22]
    charges_df = qa.predict.final_charge_dataset("all_charges.xls", "template.pdb", mutations)
