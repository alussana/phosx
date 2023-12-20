#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import h5py
from random import shuffle
from multiprocessing import Pool

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

    """
    # transform pssm into probabilities, score sequence by multiplying values and scale by prob of random peptide
    pssm_df = pssm_df.apply(lambda x: x / x.sum(), axis=1)
    n_aa = len(AA_LIST)
    n_pos = len(seq_str)
    n_res = 0
    p = 1
    for i in range(n_pos):
        if seq_str[i] != '_':
            n_res = n_res + 1
            pos = list(pssm_df.index)[i]
            p = p * pssm_df.loc[pos, seq_str[i]] * 20
    """

    """
    # add scores instead of multiplying
    p = 0
    if pssm_df.loc[0, seq_str[5]] != 0:
        for i in range(len(seq_str)):
            if seq_str[i] != '_':
                pos = list(pssm_df.index)[i]
                p = p + pssm_df.loc[pos, seq_str[i]]
    """

    return p


def pssm_scoring(seq: str, pssm_df_dict: dict):
    record = {}

    for kinase in pssm_df_dict.keys():
        p = score_sequence(seq, pssm_df_dict[kinase])
        record[kinase] = p

    out_series = pd.Series(record)

    return out_series


def binarise_pssm_scores(scaled_scores: pd.DataFrame, n: int = 15):
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


# take the largest deviation between test and null curves
def ks_statistic(
    test_var_series: pd.Series,
    null_var_df: pd.DataFrame,
    plot_bool: bool = False,
    out_plot_dir_str: str = "phosx_output",
):
    kinase = test_var_series.name

    deltas = []
    running_sum_test = 0
    running_sum_null = 0

    for i in range(len(test_var_series)):
        running_sum_test = running_sum_test + test_var_series[i]
        running_sum_null = running_sum_null + null_var_df[kinase][i]
        deltas.append(running_sum_test - running_sum_null)

    max_delta = max(deltas)
    min_delta = min(deltas)

    if plot_bool:
        data = pd.Series(deltas)

        import matplotlib.pyplot as plt
        import seaborn as sns

        if not os.path.exists(out_plot_dir_str):
            os.makedirs(out_plot_dir_str)

        plt.clf()
        fig, ax = plt.subplots(figsize=(4, 3))

        ax = sns.lineplot(data=data, linewidth=2.5)

        ax.set(
            title=f"{kinase}", xlabel="Phosphosite rank", ylabel="Running Sum Statistic"
        )

        sns.despine()
        plt.tight_layout()

        plt.savefig(os.path.join(out_plot_dir_str, f"{kinase}.svg"))
        plt.close()

    if abs(max_delta) > abs(min_delta):
        return max_delta

    else:
        return min_delta


"""
# take the deviation at the end of the test and null curves
def ks_statistic(
    
    test_var_series:pd.Series,
    null_var_df:pd.DataFrame,
    plot_bool:bool=False,
    out_plot_dir_str:str=''

):
    
    kinase = test_var_series.name
    
    running_sum_test = 0
    running_sum_null = 0
    
    running_sum_test_list = [0]
    running_sum_null_list = [0]
    
    for i in range(len(test_var_series)):
        
        running_sum_test = running_sum_test + test_var_series[i]
        running_sum_null = running_sum_null + null_var_df[kinase][i]
        
        running_sum_test_list.append(running_sum_test)
        running_sum_null_list.append(running_sum_null)
   
    if plot_bool:
        
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        if not os.path.exists(out_plot_dir_str):
            os.makedirs(out_plot_dir_str)
        
        data = pd.DataFrame({
            f'{kinase}': running_sum_test_list,
            'Null': running_sum_null_list
        })
        
        plt.clf()
        fig, ax = plt.subplots(figsize=(4, 3))
        
        ax = sns.lineplot(
            data=pd.DataFrame(data),
            linewidth=2.5
        )
        
        ax.set(
            title=f'{kinase}',
            xlabel='Phosphosite rank',
            ylabel='Statistic'
        )
        
        sns.despine()
        plt.tight_layout()
        
        plt.savefig(os.path.join(out_plot_dir_str, f'{kinase}.svg'))
        plt.close()
    
    delta = running_sum_test - running_sum_null
    
    return delta
"""


def compute_null_ks(
    y: pd.DataFrame,
    y_null: pd.DataFrame,
):
    # permute phosphosites
    idx_list = list(y.index)
    shuffle(idx_list)
    shuffled_y = y.reindex(idx_list)
    shuffled_y.index = list(range(len(shuffled_y)))

    # compute ks between shuffled pssm scores and null scores
    ks_series = shuffled_y.apply(ks_statistic, args=[y_null, False], axis=0)

    ks_df = pd.DataFrame([ks_series])

    return ks_df


def uniform_binary_vector(v: pd.Series):
    n = v.sum()

    if n == 0:
        return v

    l = len(v)
    flipped = False

    if n > (l / 2):
        v = v.apply(lambda bit: ~bit & 1)
        l = len(v)
        n = v.sum()
        flipped = True

    i_spacing = int(l / n)
    modulo = l % i_spacing

    i_start = int((modulo + i_spacing) / 2)

    uniform_v = [0 for i in range(i_start)]
    uniform_v.append(1)

    counter = n - 1
    spacing = 0
    for i in range(i_start + 1, l):
        spacing = spacing + 1
        if spacing == i_spacing and counter != 0:
            uniform_v.append(1)
            spacing = 0
            counter = counter - 1
        else:
            uniform_v.append(0)

    uniform_v = pd.Series(uniform_v)
    uniform_v.name = v.name

    if flipped:
        uniform_v = uniform_v.apply(lambda bit: ~bit & 1)

    return uniform_v


def compute_ks_empirical_distrib(
    binarised_pssm_scores: pd.DataFrame,
    seqrnk_series: pd.Series,
    n: int = 1000,
    n_proc: int = 1,
):
    # uniformly assign the binarised pssm scores along the ranks
    null_binarised_pssm_scores = binarised_pssm_scores.apply(
        uniform_binary_vector, axis=0
    )

    # multiply the binarised pssm scores by the phosphosites ranking scores
    y = binarised_pssm_scores.apply(lambda x: x * seqrnk_series, axis=0)
    y_null = null_binarised_pssm_scores.apply(lambda x: x * seqrnk_series, axis=0)

    arg1 = [y for i in range(n)]
    arg2 = [y_null for i in range(n)]

    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(compute_null_ks, zip(arg1, arg2))

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

        num_int = sum([x > value_float for x in empirical_list])
        den_int = len(empirical_list)
        pvalue = num_int / den_int

        if pvalue > 0.5:
            num_int = den_int - num_int
            pvalue = num_int / den_int
        ks_pvalue_series[kinase_str] = pvalue

        if plot_bool:
            plt.clf()
            fig, ax = plt.subplots(figsize=(5, 3))

            ax = sns.displot(data=empirical_list, height=3, aspect=5 / 3, bins=100)

            ax.set(
                title=f"Empirical ks distribution for {kinase_str} ({len(empirical_list)} perm.)",
                xlabel="KS",
                ylabel="Count",
            )

            plt.axvline(x=value_float)

            sns.despine()
            plt.tight_layout()

            plt.savefig(os.path.join(out_plot_dir_str, f"{kinase_str}_ks_dist.svg"))
            plt.close()

    return ks_pvalue_series


def compute_activity_score(results_df: pd.DataFrame):
    # compute FPR
    results_df["FPR p value"] = results_df["p value"] * len(results_df)
    results_df["FPR p value"].loc[results_df["FPR p value"] > 1] = 1

    # compute uncorrected activity score as -log10(FPR), signed with the sign of KS
    activity_score_series = results_df["p value"].copy()
    # 0. compute maximum absolute activity score
    min_non0_fpr_float = min(activity_score_series.loc[activity_score_series != 0])
    max_abs_activity_score_float = -np.log10(min_non0_fpr_float)
    # 1. replace zeroes with tiny numbers
    activity_score_series.loc[activity_score_series == 0] = np.nextafter(
        np.float32(0), np.float32(1)
    )
    # 2. compute -log10
    activity_score_series = -np.log10(activity_score_series)
    # 3. cap at maximum
    activity_score_series.loc[
        activity_score_series > max_abs_activity_score_float
    ] = max_abs_activity_score_float
    # 4. apply the same sign as ks
    results_df["uncorrected Activity Score"] = activity_score_series * results_df[
        "KS"
    ].apply(lambda x: 1 if x > 0 else -1)
    # 5. set Activity Score to 0 if direction of enrichment is not consistent with sign of KS
    results_df["Activity Score"] = results_df.apply(
        lambda x: x["uncorrected Activity Score"]
        if x["KS"] * np.sign(x["KS"]) > x["null KS mean"] * np.sign(x["KS"])
        else 0,
        axis=1,
    )
    return results_df


def kinase_activities(
    seqrnk_file: str,
    pssm_h5_file: str,
    pssm_score_quantiles_h5_file: str,
    n_perm: int = 10000,
    n_top_kinases: int = 15,
    n_proc: int = 1,
    plot_figures: bool = False,
    out_plot_dir: str = "phosx_output",
):
    pssm_df_dict = read_pssms(pssm_h5_file)
    seqrnk = read_seqrnk(seqrnk_file)
    pssm_bg_scores_df = read_pssm_score_quantiles(pssm_score_quantiles_h5_file)

    # score phosphosite sequences with each PSSM
    seq_series = seqrnk["Sequence"]
    arg1 = [seq_series[i] for i in range(len(seq_series))]
    arg2 = [pssm_df_dict for i in range(len(seq_series))]
    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(pssm_scoring, zip(arg1, arg2))
    pssm_scoring_df = pd.concat(dfs_list, axis=1).T
    pssm_scoring_df.index = list(range(len(seq_series)))

    # quantile-scale the PSSM scores for each kinase
    pssm_scoring_scaled01_df = pssm_scoring_df.apply(
        quantile_scaling, args=[pssm_bg_scores_df], axis=0
    )

    # binarise PSSM scores
    binarised_pssm_scores = binarise_pssm_scores(
        pssm_scoring_scaled01_df, n=n_top_kinases
    )

    # compute empirical distribution of ks statistic for all kinases
    ks_empirical_distrib_df = compute_ks_empirical_distrib(
        binarised_pssm_scores=binarised_pssm_scores,
        seqrnk_series=seqrnk["Score"],
        n=n_perm,
        n_proc=n_proc,
    )
    ks_empirical_mean_series = ks_empirical_distrib_df.apply(lambda x: x.mean(), axis=0)
    ks_empirical_mean_series.name = "null KS mean"

    # compute real ks statistic for all kinases
    # 1) multiply the pssm scores by the phosphosites seqrnk
    y = binarised_pssm_scores.apply(lambda x: x * seqrnk["Score"], axis=0)
    # 2) multiply the null scores by the phosphosites seqrnk
    null_binarised_pssm_scores = binarised_pssm_scores.apply(
        uniform_binary_vector, axis=0
    )
    y_null = null_binarised_pssm_scores.apply(lambda x: x * seqrnk["Score"], axis=0)
    # 3) compute ks statistic for all kinases
    ks_series = y.apply(ks_statistic, args=[y_null, plot_figures, out_plot_dir], axis=0)
    ks_series.name = "KS"

    # compute ks pvalues
    ks_pvalue_series = compute_ks_pvalues(
        ks_empirical_distrib_df=ks_empirical_distrib_df,
        ks_series=ks_series,
        plot_bool=plot_figures,
        out_plot_dir_str=out_plot_dir,
    )

    # output table
    results_df = pd.concat(
        [ks_series, ks_empirical_mean_series, ks_pvalue_series], axis=1
    )

    # compute activity score for all kinases
    results_df = compute_activity_score(results_df)

    # export results
    print(results_df.to_csv(sep="\t", header=True, index=True))
