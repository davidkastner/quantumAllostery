"""Backpropogated feature predictions from ML models."""

import pandas as pd
import numpy as np
from typing import List, Tuple
from demystifying import feature_extraction as fe, visualization
from demystifying import relevance_propagation as relprop
from sklearn.preprocessing import MinMaxScaler
from sklearn.utils import shuffle
import qa.process

def shuffle_data(samples):
    # Shuffle the data in groups of 100
    print("   > Shuffling data.")
    n_samples = samples.shape[0]
    n_samples = int(n_samples / 100) * 100
    inds = np.arange(n_samples)
    inds = inds.reshape((int(n_samples / 100), 100))
    perm_inds = np.random.permutation(inds)
    perm_inds = np.ravel(perm_inds)
    samples = samples[perm_inds]  # Apply the shuffling to the matrix
    labels = labels[perm_inds]  # Apply the same shuffling to the labels

    return samples, labels
    

def create_combined_csv(charge_files: List[str], templates: List[str], mutations: List[int]) -> pd.DataFrame:
    """
    Generate a pd.DataFrame of all features.

    Returns
    -------
    charges_df: pd.DataFrame
        The original charge data as a pandas dataframe.
    samples_df: pd.DataFrame
        One-hot-encoded labels for each frame.
    """
    
    # Convert the input data files to pd.DataFrames
    dataframes = []
    labels_df = pd.DataFrame()
    label_index = 0
    for charge_file,template in zip(charge_files,templates):
        print(f"   > Converting atoms to residues for {charge_file}.")
        # Average the charges by residue
        # We does this to minimize the inaccuracies of mulliken charges
        avg_by_residues = qa.process.average_by_residues(charge_file, template)

        # Add a column for the labels for each frame
        # You might think it would be better to use one-hot-encoding
        # However, this package specifically asks for indices for each group
        print(f"   > Creating labels for {charge_file}.")
        label = [label_index for x in range(len(avg_by_residues))]
        avg_by_residues["Label"] = label

        # Store the labeled, averaged frames in the list
        dataframes.append(avg_by_residues)

        label_index += 1 # Give the next group a different label

    # Drop the residue columns that were mutated
    # We can't compare these residues' charges as their atom counts differ
    clean_dataframes = []
    for index,df in enumerate(dataframes):
        df = df.drop(df.columns[mutations], axis=1)
        clean_dataframes.append(df)

    # Combine both dataframes into a single dataframe
    print(f"   > Combining {charge_files}.")
    combined_df = pd.concat(clean_dataframes, ignore_index=True, sort=False)
    combined_df.fillna(0, inplace=True)

    # Break off the last columns (labels) and save them as their own df
    charges_df = combined_df.iloc[:,:-1]
    labels_df = combined_df.iloc[:,-1:]

    return charges_df, labels_df


def data_processing(df, samples):
    """
    Scales the data for the ML workflows.

    Parameters
    ----------
    df: pd.DataFrame
        The original data as a pandas dataframe

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
    samples = samples.values.tolist()

    return data_norm, samples


def run_ml(data_norm, samples):
    """
    ML analysis workflow.

    Parameters
    ----------
    df_norm: numpy matrix
        The data scaled by column.
    samples_df: pd.DataFrame
        One-hot-encoded labels for each frame.

    """
    print("   > Creating csv's to check.")
    np.savetxt("data.csv", data_norm, delimiter=",")
    np.savetxt("samples.csv", samples, fmt='%i', delimiter=",")


    # Set the arguments for the ML workflows
    kwargs = {
        "samples": data_norm,
        "labels": samples,
        "filter_by_distance_cutoff": False,
        "lower_bound_distance_cutoff": 1.0,
        "upper_bound_distance_cutoff": 1.0,
        "use_inverse_distances": False,
        "n_splits": 3,
        "n_iterations": 5,
        "scaling": True,
    }

    # Running various ML workflows
    models = ["RF", "KL", "MLP"]
    feature_extractors = [
        fe.RandomForestFeatureExtractor(
            one_vs_rest=True, classifier_kwargs={"n_estimators": 100}, **kwargs
        ),
        fe.KLFeatureExtractor(**kwargs),
        fe.MlpFeatureExtractor(
            classifier_kwargs={
                "hidden_layer_sizes": (120,),
                "solver": "adam",
                "max_iter": 1000000,
            },
            activation=relprop.relu,
            **kwargs,
        ),
    ]

    # Process the results
    postprocessors = []
    working_dir = "."
    for extractor, model in zip(feature_extractors, models):
        print(f"   > Running {model} model.")
        extractor.extract_features()
        # Post-process data (rescale and filter feature importances)
        postprocessor = extractor.postprocessing(
            working_dir=working_dir,
            rescale_results=True,
            filter_results=False,
            feature_to_resids=None,
        )
        postprocessor.average()
        postprocessor.persist()
        postprocessors.append(postprocessor)

    # Create a summarizing plot of the results of the ML models
    visualization.visualize(
        [postprocessors],
        show_importance=True,
        show_projected_data=False,
        show_performance=False,
        highlighted_residues=[23],
        outfile="./importance.pdf",
    )

# Execute the script
if __name__ == "__main__":
    create_combined_csv(["mc6.xls","mc6s.xls","mc6sa.xls"], ["mc6.pdb","mc6s.pdb","mc6sa.pdb"], [0,2,15,16,19,22,27])