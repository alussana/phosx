#!/usr/bin/env python3

import os
import pandas as pd
from multiprocessing import Pool
from tqdm import tqdm


def ks_statistic(
    deltas_series: pd.Series,
    plot_bool: bool = False,
    out_plot_dir_str: str = "phosx_output",
):
    kinase = deltas_series.name

    cumsum = deltas_series.cumsum()
    max_ks = cumsum.max()
    min_ks = cumsum.min()

    if plot_bool:
        data = pd.Series(cumsum)

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
        plt.close(fig)

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
    running_sum_deltas_df = P_hit_df.copy()
    for kinase in running_sum_deltas_df.columns:
        running_sum_deltas_df[kinase].loc[running_sum_deltas_df[kinase] == 0] = (
            P_miss_series[kinase]
        )

    # compute ks for each kinase
    ks_series = running_sum_deltas_df.apply(
        ks_statistic, args=[plot_bool, out_plot_dir_str], axis=0
    )

    return ks_series


def compute_null_ks(seqrnk_series: pd.Series, binarised_pssm_scores: pd.DataFrame):
    # Shuffle the rows of binarised_pssm_scores while keeping the index aligned
    shuffled_binarised_pssm_scores = binarised_pssm_scores.sample(frac=1).reset_index(drop=True)
    
    # Compute number of non-hits for each kinase
    Nnh_series = (shuffled_binarised_pssm_scores == 0).sum(axis=0)

    seqrnk_abs_series = seqrnk_series.abs()
    
    # Compute ranking metric for hits
    rj_df = shuffled_binarised_pssm_scores.mul(seqrnk_abs_series, axis=0)
    
    # Scale ranking metric for hits to sum to 1 for each kinase
    P_hit_df = rj_df / rj_df.sum(axis=0)
    
    # Compute decrement score for non-hits
    P_miss_series = -1 / Nnh_series
    
    # Compute running sum deltas: P_hit where hit, P_miss where miss
    running_sum_deltas_df = P_hit_df.where(P_hit_df != 0, P_miss_series, axis=1)
    
    # Compute KS statistic for each kinase
    ks_series = running_sum_deltas_df.apply(ks_statistic, axis=0)
    
    return pd.DataFrame([ks_series])


def compute_ks_empirical_distrib(
    binarised_pssm_scores: pd.DataFrame,
    seqrnk_series: pd.Series,
    n: int = 1000,
    n_proc: int = 1,
):
    # run compute_null_ks n times in parallel
    arg1 = [seqrnk_series for i in range(n)]
    arg2 = [binarised_pssm_scores for i in range(n)]

    # with Pool(processes=n_proc) as pool:
    #    dfs_list = pool.starmap(compute_null_ks, zip(arg1, arg2))

    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(
            compute_null_ks,
            tqdm(
                zip(arg1, arg2), total=n, desc="    Performing permutations   ", ncols=80
            ),
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
            plt.close(fig)

    return ks_pvalue_series
