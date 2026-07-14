#!/usr/bin/env python3

import numpy as np
import pandas as pd
from phosx.cli import phosx
from phosx.kinase_activity import compute_kinase_activities
from phosx.pssm_enrichment import compute_ks_pvalues
from os import path


def test_kinase_activities_1core():
    phosx(
        seqrnk_file=str(
            path.join(
                path.dirname(__file__), "seqrnk/koksal2018_log2.fold.change.8min.seqrnk"
            )
        ),
        s_t_pssm_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/S_T_PSSMs.h5")
        ),
        y_pssm_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/Y_PSSMs.h5")
        ),
        s_t_pssm_score_quantiles_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/S_T_PSSM_score_quantiles.h5")
        ),
        y_pssm_score_quantiles_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/Y_PSSM_score_quantiles.h5")
        ),
        n_perm=10000,
        n_proc=1,
        plot_figures=False,
        out_path="phosx_output/out_kinase_activities_1core.tsv",
    )


def test_kinase_activities_4cores():
    phosx(
        seqrnk_file=str(
            path.join(
                path.dirname(__file__), "seqrnk/koksal2018_log2.fold.change.8min.seqrnk"
            )
        ),
        s_t_pssm_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/S_T_PSSMs.h5")
        ),
        y_pssm_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/Y_PSSMs.h5")
        ),
        s_t_pssm_score_quantiles_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/S_T_PSSM_score_quantiles.h5")
        ),
        y_pssm_score_quantiles_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/Y_PSSM_score_quantiles.h5")
        ),
        n_perm=10000,
        n_proc=4,
        plot_figures=False,
        out_path="phosx_output/out_kinase_activities_4cores.tsv",
    )


def test_kinase_activities_w_figures():
    phosx(
        seqrnk_file=str(
            path.join(
                path.dirname(__file__), "seqrnk/koksal2018_log2.fold.change.8min.seqrnk"
            )
        ),
        s_t_pssm_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/S_T_PSSMs.h5")
        ),
        y_pssm_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/Y_PSSMs.h5")
        ),
        s_t_pssm_score_quantiles_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/S_T_PSSM_score_quantiles.h5")
        ),
        y_pssm_score_quantiles_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/Y_PSSM_score_quantiles.h5")
        ),
        n_perm=10000,
        n_proc=4,
        plot_figures=True,
        out_path="phosx_output/out_kinase_activities_w_figures.tsv",
    )


def test_too_few_phosphosites_to_infer_activity():
    # Regression test: when the input has too few phosphosites of a given type,
    # no kinase reaches the minimum number of associated phosphosites and none
    # can be tested. Every kinase should be reported with undefined activity
    # and a valid substrates DataFrame should be returned.
    results_df, substrates_df = compute_kinase_activities(
        seqrnk_file=str(
            path.join(path.dirname(__file__), "seqrnk/few_tyr_phosphosites.seqrnk")
        ),
        pssm_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/Y_PSSMs.h5")
        ),
        pssm_score_quantiles_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/Y_PSSM_score_quantiles.h5")
        ),
        n_perm=100,
        n_proc=1,
        tyr_only=True,
    )

    # nothing could be inferred, but all kinases are reported: KS and p values
    # are undefined and the activity defaults to 0 (as for any untested kinase)
    assert len(results_df) > 0
    assert results_df["KS"].isna().all()
    assert results_df["p value"].isna().all()
    assert (results_df["Activity Score"] == 0).all()

    # substrates must be a DataFrame (not None) so downstream steps keep working
    assert substrates_df is not None
    assert list(substrates_df.columns)


def test_pvalue_never_zero_and_capped_at_resolution():
    # Regression test: a permutation p value must never be reported as 0. With
    # the (b + 1) / (m + 1) estimator, when the observed KS statistic is more
    # extreme than all m permutations the p value hits its floor 1 / (m + 1),
    # not 0. Cover both the positive- and negative-KS branches.
    n_perm = 100
    empirical_distrib_df = pd.DataFrame(
        {"KINASE1": np.linspace(-0.5, 0.5, n_perm)}
    )

    # observed KS larger than every permutation -> minimum attainable p value
    ks_series = pd.Series({"KINASE1": 10.0})
    ks_series.name = "KS"
    pvals = compute_ks_pvalues(empirical_distrib_df, ks_series)
    assert (pvals > 0).all()
    assert (pvals <= 1).all()
    assert pvals["KINASE1"] == 1 / (n_perm + 1)

    # observed KS smaller than every permutation -> same floor via the < branch
    ks_series = pd.Series({"KINASE1": -10.0})
    ks_series.name = "KS"
    pvals = compute_ks_pvalues(empirical_distrib_df, ks_series)
    assert (pvals > 0).all()
    assert pvals["KINASE1"] == 1 / (n_perm + 1)


def test_pvalues_are_bounded_in_full_run():
    # Regression test on a full run: every reported p value must be strictly
    # positive, at most 1, and at least the resolution floor 1 / (n_perm + 1),
    # so that no zero p value is ever written to the output.
    n_perm = 100
    results_df, _ = compute_kinase_activities(
        seqrnk_file=str(
            path.join(
                path.dirname(__file__), "seqrnk/koksal2018_log2.fold.change.8min.seqrnk"
            )
        ),
        pssm_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/S_T_PSSMs.h5")
        ),
        pssm_score_quantiles_h5_file=str(
            path.join(path.dirname(__file__), "../phosx/data/S_T_PSSM_score_quantiles.h5")
        ),
        n_perm=n_perm,
        n_proc=1,
        ser_thr_only=True,
    )

    p = results_df["p value"].dropna()
    assert len(p) > 0
    assert (p > 0).all()
    assert (p <= 1).all()
    assert (p >= 1 / (n_perm + 1)).all()


def test_all():
    test_kinase_activities_1core()
    test_kinase_activities_4cores()
    test_kinase_activities_w_figures()
    test_too_few_phosphosites_to_infer_activity()
    test_pvalue_never_zero_and_capped_at_resolution()
    test_pvalues_are_bounded_in_full_run()


if __name__ == "__main__":
    test_all()
