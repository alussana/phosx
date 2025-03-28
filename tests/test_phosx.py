#!/usr/bin/env python3

from phosx.cli import phosx
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


def test_all():
    test_kinase_activities_1core()
    test_kinase_activities_4cores()
    test_kinase_activities_w_figures()


if __name__ == "__main__":
    test_all()
