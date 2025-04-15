#!/usr/bin/env python3

import pandas as pd
import numpy as np


def quantile_scaling(x: pd.Series, sorted_bg_scores_dict: dict):
    """
    Compute quantile scores for the given series x using precomputed sorted background scores.

    Args:
        x (pd.Series): Series of scores for a specific kinase.
        sorted_bg_scores_dict (dict): Dictionary with kinase names as keys and sorted numpy arrays of background scores as values.

    Returns:
        pd.Series: Series of quantile scores corresponding to x.
    """
    sorted_bg_scores = sorted_bg_scores_dict[x.name]
    x_values = x.values

    num = np.searchsorted(sorted_bg_scores, x_values, side="right")
    den = len(sorted_bg_scores)
    quantile_scores = num / den

    scores = pd.Series(quantile_scores, index=x.index)
    scores.name = x.name

    return scores


def score_sequence(seq_str: str, pssm_df: pd.DataFrame):
    if len(seq_str) != len(pssm_df):
        raise Exception("Sequence length cannot be different from pssm length")

    # score sequence by multiplying values and scale by prob of random peptide as described in Supplementary Note 2 of https://doi.org/10.1038/s41586-022-05575-3
    n_pos = len(seq_str)
    p = 1
    for i in range(n_pos):
        if seq_str[i] != "_":
            try:
                pos = list(pssm_df.index)[i]
                p = p * pssm_df.loc[pos, seq_str[i]]
            except KeyError:
                print("Non-canonical amino acid symbol found in sequence.")
                print(
                    f"'{seq_str[i]}' was found, but kinase PSSMs can handle the following symbols:"
                )
                print(f"{' '.join(AA_LIST)}")
    return p


def max_score_sequence(seq_str: str, pssm_df: pd.DataFrame):
    """
    Compute the highest PSSM score among all possible substrings of length equal to pssm_df length.

    Args:
        seq_str (str): Input sequence string of length >= len(pssm_df)
        pssm_df (pd.DataFrame): Position-specific scoring matrix DataFrame

    Returns:
        float: Highest score among all possible substrings
    """
    pssm_length = len(pssm_df)
    if len(seq_str) < pssm_length:
        raise Exception("Sequence length cannot be shorter than pssm length")
    max_score = float("-inf")

    # Iterate through all possible substrings of length pssm_length
    for i in range(len(seq_str) - pssm_length + 1):
        substring = seq_str[i : i + pssm_length]
        p = score_sequence(substring, pssm_df)
        max_score = max(max_score, p)

    return max_score


def pssm_scoring(seq: str, pssm_df_dict: dict):
    record = {}

    for kinase in pssm_df_dict.keys():
        p = score_sequence(seq, pssm_df_dict[kinase])
        record[kinase] = p

    out_series = pd.Series(record)

    return out_series


def max_pssm_scoring(seq: str, pssm_df_dict: dict):
    record = {}

    for kinase in pssm_df_dict.keys():
        p = max_score_sequence(seq, pssm_df_dict[kinase])
        record[kinase] = p

    out_series = pd.Series(record)

    return out_series


def binarise_pssm_scores(scaled_scores: pd.DataFrame, n: int = 5, m: float = 0.95):
    """Binarise kinase PSSM scores given the number of top-scoring kinases that should be assigned to a phosphosite and the minimum PSSM score quantile.

    Args:
        scaled_scores (pandas.DataFrame): kinase PSSM scores for a list of phosphosites; rows are phosphosites, columns are kinases.
        n (int, optional): number of top scoring kinases to assign to each phopshosite. Defaults to 5.
        m (float, optional): minimum PSSM score quantile that a phosphosite has to satisfy to be potentially assigned to a kinase. Defaults to 0.95.

    Returns:
        pandas.DataFrame: binarised kinase PSSM scores for the given phosphosites.
    """
    scaled_scores = scaled_scores.copy()

    def find_bin_threshold(series: pd.Series):

        sorted_values = sorted(series.values, reverse=True)
        threshold = sorted_values[n]
        if threshold < m:
            threshold = m
        series.loc[series > threshold] = 1
        series.loc[series <= threshold] = 0

        return series

    binarised_scores = scaled_scores.apply(find_bin_threshold, axis=1)

    return binarised_scores
