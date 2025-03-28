#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import h5py
import math

AA_LIST = [
    "G",
    "P",
    "A",
    "V",
    "L",
    "I",
    "M",
    "C",
    "F",
    "Y",
    "W",
    "H",
    "K",
    "R",
    "Q",
    "N",
    "E",
    "D",
    "S",
    "T",
    "s",
    "t",
    "y",
]

POSITIONS_LIST = list(range(-5, 5))


def softmax(x):
    # Subtract max for numerical stability
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()


def decay_from_1(x, decay_factor=64):
    return math.exp(-(decay_factor * ((x - 1) ** 2)))


def read_pssms(pssms_h5_file: str):
    pssms_h5 = h5py.File(pssms_h5_file, "r")
    pssm_df_dict = {}

    for kinase in pssms_h5.keys():
        pssm_df_dict[kinase] = pd.DataFrame(pssms_h5[kinase])
        pssm_df_dict[kinase].columns = AA_LIST
        pssm_df_dict[kinase].index = POSITIONS_LIST

    return pssm_df_dict


def read_seqrnk(seqrnk_file: str, ser_thr_only: bool = False, tyr_only: bool = False):
    if ser_thr_only and tyr_only:
        print(
            "W: both ser_thr_only and tyr_only have been set to True. PhosX will consider only Ser/Thr phosphosites in this run.",
            file=sys.stderr,
        )

    seqrnk = pd.read_csv(seqrnk_file, sep="\t", header=None)
    seqrnk.columns = ["Sequence", "Score"]

    # remove phosphosites with non existent ranking metric
    seqrnk = seqrnk.loc[seqrnk["Score"] != 0,]
    seqrnk = seqrnk.dropna(axis=0)
    seqrnk.index = range(len(seqrnk))

    if ser_thr_only:
        seqrnk = seqrnk.loc[
            [
                seqrnk["Sequence"][i][5] == "S" or seqrnk["Sequence"][i][5] == "T"
                for i in range(len(seqrnk))
            ]
        ]
        seqrnk.index = range(len(seqrnk))
    elif tyr_only:
        seqrnk = seqrnk.loc[
            [seqrnk["Sequence"][i][5] == "Y" for i in range(len(seqrnk))]
        ]
        seqrnk.index = range(len(seqrnk))

    return seqrnk


def read_pssm_score_quantiles(pssm_score_quantiles_h5_file: str):
    pssm_bg_scores_df = pd.read_hdf(pssm_score_quantiles_h5_file, key="pssm_scores")

    return pssm_bg_scores_df


def hdf5_to_dict(hdf5_file_path):
    """
    Convert an HDF5 file's contents into a nested dictionary recursively,
    with special handling for string datasets.
    
    Args:
        hdf5_file_path (str): Path to the HDF5 file
        
    Returns:
        dict: Nested dictionary containing the HDF5 file's structure and data
    """
    def _convert_group(h5_group):
        result = {}
        for key, item in h5_group.items():
            # If item is a group, recursively convert it
            if isinstance(item, h5py.Group):
                result[key] = _convert_group(item)
            # If item is a dataset
            elif isinstance(item, h5py.Dataset):
                # Check if the dataset contains string type
                if h5py.check_string_dtype(item.dtype):
                    # Handle scalar string
                    if item.shape == ():
                        value = item[()]
                        result[key] = value.decode('utf-8') if isinstance(value, bytes) else str(value)
                    # Handle array of strings
                    else:
                        value = item[:]
                        result[key] = [v.decode('utf-8') if isinstance(v, bytes) else str(v) 
                                     for v in value]
                else:
                    # Handle non-string scalar datasets
                    if item.shape == ():
                        result[key] = item[()]
                    # Handle non-string array datasets
                    else:
                        result[key] = item[:]
                    # Convert numpy types to Python native types if possible
                    if hasattr(result[key], 'tolist'):
                        result[key] = result[key].tolist()
        return result
    
    # Open the HDF5 file and convert its contents
    try:
        with h5py.File(hdf5_file_path, 'r') as f:
            return _convert_group(f)
    except Exception as e:
        raise Exception(f"Error reading HDF5 file: {str(e)}")
    

def sync_dataframe_with_names(df, names_list, fill_value=0):
    """
    Synchronize a DataFrame's rows and columns with a list of names.
    
    Parameters:
    - df (pd.DataFrame): Input DataFrame
    - names_list (list): List of strings to sync rows and columns with
    - fill_value (scalar, optional): Value to fill for new entries (default: 0)
    
    Returns:
    - pd.DataFrame: Synchronized DataFrame
    """
    # Convert names_list to set for comparison
    names_set = set(names_list)
    current_cols = set(df.columns)
    current_idx = set(df.index)
    
    # Create a new DataFrame with names_list as both index and columns
    new_df = pd.DataFrame(fill_value, 
                         index=names_list, 
                         columns=names_list)
    
    # Find common names between df and names_list for both rows and columns
    common_rows = current_idx & names_set
    common_cols = current_cols & names_set
    
    # Copy data from the original DataFrame for overlapping rows and columns
    if common_rows and common_cols:
        new_df.loc[list(common_rows), list(common_cols)] = df.loc[list(common_rows), list(common_cols)]
    
    # Ensure the order matches names_list for both rows and columns
    new_df = new_df.loc[names_list, names_list]
    
    return new_df

