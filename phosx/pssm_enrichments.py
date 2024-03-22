#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import h5py
from random import shuffle
from multiprocessing import Pool
from tqdm import tqdm

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


def read_pssms(pssms_h5_file: str):
    pssms_h5 = h5py.File(pssms_h5_file, "r")
    pssm_df_dict = {}

    for kinase in pssms_h5.keys():
        pssm_df_dict[kinase] = pd.DataFrame(pssms_h5[kinase])
        pssm_df_dict[kinase].columns = AA_LIST
        pssm_df_dict[kinase].index = POSITIONS_LIST

    return pssm_df_dict


def read_seqrnk(seqrnk_file: str):
    seqrnk = pd.read_csv(seqrnk_file, sep="\t", header=None)
    seqrnk.columns = ["Sequence", "Score"]

    # remove phosphosites with non existent ranking metric
    seqrnk = seqrnk.loc[seqrnk["Score"] != 0,]
    seqrnk = seqrnk.dropna(axis=0)
    seqrnk.index = range(len(seqrnk))

    # tmp: only consider S/T phosphosites
    seqrnk = seqrnk.loc[
        [
            seqrnk["Sequence"][i][5] == "S" or seqrnk["Sequence"][i][5] == "T"
            for i in range(len(seqrnk))
        ]
    ]
    seqrnk.index = range(len(seqrnk))

    return seqrnk


def read_pssm_score_quantiles(pssm_score_quantiles_h5_file: str):
    pssm_bg_scores_df = pd.read_hdf(pssm_score_quantiles_h5_file, key="pssm_scores")

    return pssm_bg_scores_df


def quantile_scaling(x: pd.Series, pssm_bg_scores_df: pd.DataFrame):
    pssm_bg_scores_series = pssm_bg_scores_df[x.name]
    quantile_scores_list = []

    den = len(pssm_bg_scores_series)

    for i in range(len(x)):
        num = (pssm_bg_scores_series <= x[i]).sum()
        quantile_scores_list.append(num / den)

    scores = pd.Series(quantile_scores_list)
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
            pos = list(pssm_df.index)[i]
            p = p * pssm_df.loc[pos, seq_str[i]]

    return p


def pssm_scoring(seq: str, pssm_df_dict: dict):
    record = {}

    for kinase in pssm_df_dict.keys():
        p = score_sequence(seq, pssm_df_dict[kinase])
        record[kinase] = p

    out_series = pd.Series(record)

    return out_series


def binarise_pssm_scores(scaled_scores: pd.DataFrame, n: int = 8):
    """Binarise kinase PSSM scores given the number of top-scoring kinases that should be assigned to a phosphosite.

    Args:
        scaled_scores (pandas.DataFrame): kinase PSSM scores for a list of phosphosites; rows are phosphosites, columns are kinases.
        n (int, optional): number of top scoring kinases to assign to each phopshosite. Defaults to 15.

    Returns:
        pandas.DataFrame: binarised kinase PSSM scores for the given phosphosites.
    """

    def find_bin_threshold(series: pd.Series):
        sorted_values = sorted(series.values, reverse=True)
        threshold = sorted_values[n]

        series.loc[series > threshold] = 1
        series.loc[series <= threshold] = 0

        return series

    binarised_scores = scaled_scores.apply(find_bin_threshold, axis=1)

    return binarised_scores


def ks_statistic(
    deltas_series: pd.Series,
    plot_bool: bool = False,
    out_plot_dir_str: str = "phosx_output",
):
    kinase = deltas_series.name

    running_sum = 0
    running_sum_to_i = [0]

    for i in range(len(deltas_series)):
        running_sum = running_sum + deltas_series.iloc[i]
        running_sum_to_i.append(running_sum)

    max_ks = max(running_sum_to_i)
    min_ks = min(running_sum_to_i)

    if plot_bool:
        data = pd.Series(running_sum_to_i)

        import matplotlib.pyplot as plt
        import seaborn as sns

        if not os.path.exists(out_plot_dir_str):
            os.makedirs(out_plot_dir_str)

        plt.clf()
        fig, ax = plt.subplots(figsize=(4, 3))

        ax = sns.lineplot(data=data, linewidth=1.5, color="black")

        ax.set(title=f"{kinase}", xlabel="Phosphosite rank", ylabel="Running sum")

        plt.axhline(
            y=0,
            lw=0.5,
            color="black",
            linestyle="--",
        )

        sns.despine()
        plt.tight_layout()

        plt.savefig(os.path.join(out_plot_dir_str, f"{kinase}.pdf"))
        plt.close()

    if abs(max_ks) > abs(min_ks):
        return max_ks

    else:
        return min_ks


def compute_ks(
    seqrnk_series: pd.Series,
    binarised_pssm_scores: pd.DataFrame,
    plot_bool: bool = True,
    out_plot_dir_str: str = "phosx_output",
):
    # ranking metric for hits for each kinase
    rj_df = binarised_pssm_scores.apply(lambda x: x * seqrnk_series.abs(), axis=0)

    # number of non-hits for each kinase
    Nnh_series = rj_df.apply(lambda x: len(x.loc[x == 0]), axis=0)

    # scale ranking metric for hits in order to sum to 1 for each kinase
    P_hit_df = rj_df.apply(lambda x: x / x.sum(), axis=0)

    # assign decrement score for non-hits in order to sum to -1 for each kinase
    P_miss_series = -1 / Nnh_series

    # make table of running sum deltas for each kinase
    running_sum_deltas_df = P_hit_df
    for kinase in running_sum_deltas_df.columns:
        running_sum_deltas_df[kinase].loc[
            running_sum_deltas_df[kinase] == 0
        ] = P_miss_series[kinase]

    # compute ks for each kinase
    ks_series = running_sum_deltas_df.apply(
        ks_statistic, args=[plot_bool, out_plot_dir_str], axis=0
    )

    return ks_series


def compute_null_ks(seqrnk_series: pd.Series, binarised_pssm_scores: pd.DataFrame):
    # randomly permute phosphosites
    idx_list = list(seqrnk_series.index)
    shuffle(idx_list)
    shuffled_binarised_pssm_scores = binarised_pssm_scores.copy()
    shuffled_binarised_pssm_scores.index = idx_list
    shuffled_binarised_pssm_scores = shuffled_binarised_pssm_scores.sort_index()

    # ranking metric for hits for each kinase
    rj_df = shuffled_binarised_pssm_scores.apply(
        lambda x: x * seqrnk_series.abs(), axis=0
    )

    # number of non-hits for each kinase
    Nnh_series = rj_df.apply(lambda x: len(x.loc[x == 0]), axis=0)

    # scale ranking metric for hits for each kinase
    P_hit_df = rj_df.apply(lambda x: x / x.sum(), axis=0)

    # assign decrement score for non-hits in order to sum to -1 for each kinase
    P_miss_series = -1 / Nnh_series

    # make table of absolute running sum deltas for each kinase
    running_sum_deltas_df = P_hit_df
    for kinase in running_sum_deltas_df.columns:
        running_sum_deltas_df[kinase].loc[
            running_sum_deltas_df[kinase] == 0
        ] = P_miss_series[kinase]

    # compute ks for each kinase
    ks_series = running_sum_deltas_df.apply(ks_statistic, axis=0)

    ks_df = pd.DataFrame([ks_series])

    return ks_df


def compute_ks_empirical_distrib(
    binarised_pssm_scores: pd.DataFrame,
    seqrnk_series: pd.Series,
    n: int = 1000,
    n_proc: int = 1,
):
    # run compute_null_ks n times in parallel
    arg1 = [seqrnk_series for i in range(n)]
    arg2 = [binarised_pssm_scores for i in range(n)]
    
    #with Pool(processes=n_proc) as pool:
    #    dfs_list = pool.starmap(compute_null_ks, zip(arg1, arg2))
        
    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(
            compute_null_ks,
            tqdm(zip(arg1, arg2), total=n, desc='   Performing permutations   ', ncols=80)
        )

    out_df = pd.concat(dfs_list, axis=0)
    out_df.index = list(range(len(out_df)))

    return out_df


def compute_ks_pvalues(
    ks_empirical_distrib_df: pd.DataFrame,
    ks_series: pd.Series,
    plot_bool: bool = False,
    out_plot_dir_str: str = "",
):
    if plot_bool:
        import matplotlib.pyplot as plt
        import seaborn as sns

        if not os.path.exists(out_plot_dir_str):
            os.makedirs(out_plot_dir_str)

    ks_pvalue_series = pd.Series()
    ks_pvalue_series.name = "p value"

    for kinase_str in ks_empirical_distrib_df.columns:
        empirical_list = ks_empirical_distrib_df[kinase_str]

        value_float = ks_series[kinase_str]

        if value_float > 0:
            num_int = sum([x > value_float for x in empirical_list])
        else:
            num_int = sum([x < value_float for x in empirical_list])

        den_int = len(empirical_list)
        pvalue = num_int / den_int

        ks_pvalue_series[kinase_str] = pvalue

        if plot_bool:
            plt.clf()
            fig, ax = plt.subplots(figsize=(4, 3))

            ax = sns.displot(
                data=empirical_list, height=3, aspect=4 / 3, bins=100, color="white"
            )

            ax.set(
                title=f"{kinase_str} | {len(empirical_list)} permutations",
                xlabel="Enrichment score",
                ylabel="Count",
            )

            plt.axvline(x=value_float, color="black", linestyle="--")

            sns.despine()
            plt.tight_layout()

            plt.savefig(os.path.join(out_plot_dir_str, f"{kinase_str}_ks_dist.pdf"))
            plt.close()

    return ks_pvalue_series


def compute_activity_score(results_df: pd.DataFrame, max_abs_score: float):
    # compute FDR
    results_df["FDR q value"] = results_df["p value"] * len(results_df)
    results_df["FDR q value"].loc[results_df["FDR q value"] > 1] = 1

    # compute activity score as -log10(FDR), signed with the sign of KS
    activity_score_series = results_df["p value"].copy()

    # 1. replace zeroes with tiny numbers
    activity_score_series.loc[activity_score_series == 0] = np.nextafter(
        np.float32(0), np.float32(1)
    )

    # 2. compute -log10
    activity_score_series = -np.log10(activity_score_series)

    # 3. cap at maximum
    activity_score_series.loc[activity_score_series > max_abs_score] = max_abs_score

    # 4. apply the same sign as ks
    results_df["Activity Score"] = activity_score_series * np.sign(results_df["KS"])

    return results_df


def kinase_activities(
    seqrnk_file: str,
    pssm_h5_file: str,
    pssm_score_quantiles_h5_file: str,
    n_perm: int = 1000,
    n_top_kinases: int = 8,
    min_n_hits: int = 4,
    n_proc: int = 1,
    plot_figures: bool = False,
    out_plot_dir: str = "phosx_output",
    out_path = None
):
   
    print('    Loading input objects    : ', file=sys.stderr, end='')
    pssm_df_dict = read_pssms(pssm_h5_file)
    seqrnk = read_seqrnk(seqrnk_file)
    pssm_bg_scores_df = read_pssm_score_quantiles(pssm_score_quantiles_h5_file)

    # sort seqrnk by phosphosite score in descending order
    seqrnk = seqrnk.sort_values(by="Score", ascending=False)
    print('DONE', file=sys.stderr)

    # score phosphosite sequences with each PSSM
    seq_series = seqrnk["Sequence"]
    arg1 = [seq_series[i] for i in range(len(seq_series))]
    arg2 = [pssm_df_dict for i in range(len(seq_series))]
    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(
            pssm_scoring,
            tqdm(zip(arg1, arg2), ncols=80, total=len(seq_series), desc='   Scoring phosphopeptides   ')
        )
    pssm_scoring_df = pd.concat(dfs_list, axis=1).T
    pssm_scoring_df.index = list(range(len(seq_series)))

    print('Assigning kinases to targets : ', file=sys.stderr, end='')
    # quantile-scale the PSSM scores for each kinase
    pssm_scoring_scaled01_df = pssm_scoring_df.apply(
        quantile_scaling,
        args=[pssm_bg_scores_df],
        axis=0,
    )

    # binarise PSSM scores - TODO only return binarised_pssm_scores, don't change pssm_scoring_scaled01_df
    binarised_pssm_scores = binarise_pssm_scores(
        pssm_scoring_scaled01_df, n=n_top_kinases
    )

    # drop kinases with less than min_n_hits
    binarised_pssm_scores = binarised_pssm_scores.loc[
        :, binarised_pssm_scores.sum() >= min_n_hits
    ]
    print('DONE', file=sys.stderr)

    # compute empirical distribution of ks statistic for all kinases
    ks_empirical_distrib_df = compute_ks_empirical_distrib(
        binarised_pssm_scores=binarised_pssm_scores,
        seqrnk_series=seqrnk["Score"],
        n=n_perm,
        n_proc=n_proc,
    )

    # compute real ks statistic for all kinases
    ks_series = compute_ks(
        seqrnk_series=seqrnk["Score"],
        binarised_pssm_scores=binarised_pssm_scores,
        plot_bool=plot_figures,
        out_plot_dir_str=out_plot_dir,
    )
    ks_series.name = "KS"

    # compute ks p values for all kinases
    print('  Computing activity scores  : ', file=sys.stderr, end='')
    ks_pvalue_series = compute_ks_pvalues(
        ks_empirical_distrib_df=ks_empirical_distrib_df,
        ks_series=ks_series,
        plot_bool=plot_figures,
        out_plot_dir_str=out_plot_dir,
    )

    # output table
    results_df = pd.concat(
        [ks_series, ks_pvalue_series], axis=1
    )

    # compute activity score for all kinases
    results_df = compute_activity_score(results_df, np.log10(n_perm))

    # export results
    if out_path == None:
        print(results_df.to_csv(sep="\t", header=True, index=True))
    else:
        results_df.to_csv(out_path, sep="\t", header=True, index=True)
    print('DONE', file=sys.stderr, end='\n\n')
    
    return results_df
