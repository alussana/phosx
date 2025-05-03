#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm

from phosx.utils import read_pssms, read_pssm_score_quantiles, read_seqrnk
from phosx.pssms import quantile_scaling, pssm_scoring, binarise_pssm_scores
from phosx.pssm_enrichment import (
    compute_ks_empirical_distrib,
    compute_ks_pvalues,
    compute_ks,
)


def compute_activity_score(results_df: pd.DataFrame, max_abs_score: float):
    # compute FDR
    results_df["FDR q value"] = (results_df["p value"] * len(results_df)).round(
        decimals=5
    )
    # cap FDR q value at 1
    results_df.loc[results_df["FDR q value"] > 1, "FDR q value"] = 1

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


def compute_kinase_activities(
    seqrnk_file: str,
    pssm_h5_file: str,
    pssm_score_quantiles_h5_file: str,
    n_perm: int = 10000,
    n_top_kinases: int = 5,
    min_n_hits: int = 4,
    min_quantile: float = 0.95,
    n_proc: int = 1,
    plot_figures: bool = False,
    out_plot_dir: str = "phosx_output/",
    ser_thr_only=False,
    tyr_only=False,
):

    print("     Loading input objects    : ", file=sys.stderr, end="")
    pssm_df_dict = read_pssms(pssm_h5_file)
    seqrnk = read_seqrnk(seqrnk_file, ser_thr_only, tyr_only)
    pssm_bg_scores_df = read_pssm_score_quantiles(pssm_score_quantiles_h5_file)

    # sort seqrnk by phosphosite score in descending order
    seqrnk.sort_values(by="Score", ascending=False, inplace=True, ignore_index=True)
    print("DONE", file=sys.stderr)

    # score phosphosite sequences with each PSSM
    seq_series = seqrnk["Sequence"]
    arg1 = [seq_series[i] for i in range(len(seq_series))]
    arg2 = [pssm_df_dict for i in range(len(seq_series))]
    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(
            pssm_scoring,
            tqdm(
                zip(arg1, arg2),
                ncols=80,
                total=len(seq_series),
                desc="     Scoring phosphosites     ",
            ),
        )
    if len(dfs_list) == 0:
        return None
    pssm_scoring_df = pd.concat(dfs_list, axis=1).T
    pssm_scoring_df.index = list(range(len(seq_series)))

    print(" Assigning kinases to targets : ", file=sys.stderr, end="")

    # quantile-scale the PSSM scores for each kinase
    sorted_bg_scores_dict = {
        kinase: np.sort(pssm_bg_scores_df[kinase].values)
        for kinase in pssm_bg_scores_df.columns
    }
    pssm_scoring_scaled01_df = pssm_scoring_df.apply(
        quantile_scaling,
        args=[sorted_bg_scores_dict],
        axis=0,
    )

    # binarise PSSM scores
    binarised_pssm_scores = binarise_pssm_scores(
        pssm_scoring_scaled01_df, n=n_top_kinases, m=min_quantile
    )

    # drop kinases with less than min_n_hits
    binarised_pssm_scores = binarised_pssm_scores.loc[
        :, binarised_pssm_scores.sum() >= min_n_hits
    ]
    print("DONE", file=sys.stderr)

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
    ).round(decimals=5)
    ks_series.name = "KS"

    # compute ks p values for all kinases
    print("   Computing activity scores  : ", file=sys.stderr, end="")
    ks_pvalue_series = compute_ks_pvalues(
        ks_empirical_distrib_df=ks_empirical_distrib_df,
        ks_series=ks_series,
        plot_bool=plot_figures,
        out_plot_dir_str=out_plot_dir,
    ).round(decimals=5)

    # output table
    results_df = pd.concat([ks_series, ks_pvalue_series], axis=1)

    # compute activity score for all kinases
    results_df = compute_activity_score(results_df, np.log10(n_perm))
    results_df["Activity Score"] = results_df["Activity Score"].round(decimals=5)

    # add the kinases for which no inferece could be made to results_df
    kinase_list = list(pssm_df_dict.keys())
    missing_kinases = list(set(kinase_list) - set(results_df.index))
    if len(missing_kinases) > 0:
        default_values = [np.nan, np.nan, np.nan, 0]
        columns = results_df.columns
        new_rows = pd.DataFrame(
            [default_values] * len(missing_kinases),
            index=missing_kinases,
            columns=columns,
        )
        results_df = pd.concat([results_df, new_rows])

    print("DONE", file=sys.stderr, end="\n\n")

    return results_df, binarised_pssm_scores
