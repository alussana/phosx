#!/usr/bin/env python3

import argparse
from os import path, makedirs
import sys
import pandas as pd
from phosx.kinase_activity import compute_kinase_activities
from phosx.activation_evidence import compute_activation_evidence


def parse_phosx_args():
    parser = argparse.ArgumentParser(
        prog="phosx",
        description="Data-driven differential kinase activity inference from phosphosproteomics data",
        epilog="",
    )
    parser.add_argument("seqrnk", type=str, help="Path to the seqrnk file.")
    parser.add_argument(
        "-yp",
        "--y-pssm",
        type=str,
        default=str(path.join(path.dirname(__file__), "data/Y_PSSMs.h5")),
        help="Path to the h5 file storing custom Tyr PSSMs; defaults to built-in PSSMs",
    )
    parser.add_argument(
        "-stp",
        "--s-t-pssm",
        type=str,
        default=str(path.join(path.dirname(__file__), "data/S_T_PSSMs.h5")),
        help="Path to the h5 file storing custom Ser/Thr PSSMs; defaults to built-in PSSMs",
    )
    parser.add_argument(
        "-yq",
        "--y-pssm-quantiles",
        type=str,
        default=str(
            path.join(path.dirname(__file__), "data/Y_PSSM_score_quantiles.h5")
        ),
        help="Path to the h5 file storing custom Tyr kinases PSSM score quantile distributions under the key 'pssm_scores'; defaults to built-in PSSM scores quantiles",
    )
    parser.add_argument(
        "-stq",
        "--s-t-pssm-quantiles",
        type=str,
        default=str(
            path.join(path.dirname(__file__), "data/S_T_PSSM_score_quantiles.h5")
        ),
        help="Path to the h5 file storing custom Ser/Thr kinases PSSM score quantile distributions under the key 'pssm_scores'; defaults to built-in PSSM scores quantiles",
    )
    parser.add_argument(
        "-no-uae",
        "--no-upstream-activation-evidence",
        action="store_true",
        default=False,
        help="Do not compute upstream activation evidence to modify the activity scores of kinases with correlated activity; default: False",
    )
    parser.add_argument(
        "-meta",
        "--kinase-metadata",
        type=str,
        default=str(
            path.join(path.dirname(__file__), "data/kinase_metadata_annotated.h5")
        ),
        help='Path to the h5 file storing kinase metadata ("aloop_seq"); defaults to built-in metadata',
    )
    parser.add_argument(
        "-n",
        "--n-permutations",
        type=int,
        default=10000,
        help="Number of random permutations; default: 10000",
    )
    parser.add_argument(
        "-stk",
        "--s-t-n-top-kinases",
        type=int,
        default=5,
        help="Number of top-scoring Ser/Thr kinases potentially associatiated to a given phosphosite; default: 10",
    )
    parser.add_argument(
        "-yk",
        "--y-n-top-kinases",
        type=int,
        default=10,
        help="Number of top-scoring Tyr kinases potentially associatiated to a given phosphosite; default: 5",
    )
    parser.add_argument(
        "-astqth",
        "--a-loop-s-t-quantile-threshold",
        type=float,
        default=0.95,
        help="Minimum Ser/Thr PSSM score quantile for an activation loop to be considered potential substrate of a kinase; default: 0.95",
    )
    parser.add_argument(
        "-ayqth",
        "--a-loop-y-quantile-threshold",
        type=float,
        default=0.95,
        help="Minimum Tyr PSSM score quantile for an activation loop to be considered potential substrate of a kinase; default: 0.95",
    )
    parser.add_argument(
        "-urt",
        "--upreg-redundancy-threshold",
        type=float,
        default=0.5,
        help="Minimum Jaccard index of target substrates to consider two upregulated kinases having potentially correlated activity; upstream activation evidence is used to prioritize the activity of individual ones; default: 0.5",
    )
    parser.add_argument(
        "-drt",
        "--downreg-redundancy-threshold",
        type=float,
        default=0.5,
        help="Minimum Jaccard index of target substrates to consider two downregulated kinases having potentially correlated activity; upstream activation evidence is used to prioritize the activity of individual ones; default: 0.5",
    )
    parser.add_argument(
        "-mh",
        "--min-n-hits",
        type=int,
        default=4,
        help="Minimum number of phosphosites associated with a kinase for the kinase to be considered in the analysis; default: 4",
    )
    parser.add_argument(
        "-stmq",
        "--s-t-min-quantile",
        type=float,
        default=0.95,
        help="Minimum PSSM score quantile that a phosphosite has to satisfy to be potentially assigned to a Ser/Thr kinase; default: 0.95",
    )
    parser.add_argument(
        "-ymq",
        "--y-min-quantile",
        type=float,
        default=0.90,
        help="Minimum PSSM score quantile that a phosphosite has to satisfy to be potentially assigned to a Tyr kinase; default: 0.90",
    )
    parser.add_argument(
        "-df1",
        "--decay-factor",
        type=float,
        default=64,
        help="Decay factor for the exponential decay of the activation evidence when competing kinases have different activation scores. See utils.decay_from_1(); default: 64",
    )
    parser.add_argument(
        "-c",
        "--n-proc",
        type=int,
        default=1,
        help="Number of cores used for multithreading; default: 1",
    )
    parser.add_argument(
        "--plot-figures",
        action="store_true",
        help="Save figures in pdf format; see also --output_dir",
    )
    parser.add_argument(
        "-d",
        "--output-dir",
        type=str,
        default="phosx_output/",
        help="Output files directory; only relevant if used with --plot_figures; defaults to 'phosx_output/'",
    )
    parser.add_argument(
        "-nd",
        "--network-path",
        type=str,
        default=None,
        help="Output file path for the inferred Kinase-->A-loop network; if not specified, the network will not be saved",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        type=str,
        default=None,
        help="Main output table with differential kinase activitiy scores; if not specified it will be printed in STDOUT",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="0.16.0",
        help="Print package version and exit",
    )
    args = parser.parse_args()
    return args


def phosx(
    seqrnk_file: str,
    s_t_pssm_h5_file: str = str(
        path.join(path.dirname(__file__), "../phosx/data/S_T_PSSMs.h5")
    ),
    s_t_pssm_score_quantiles_h5_file: str = str(
        path.join(path.dirname(__file__), "../phosx/data/S_T_PSSM_score_quantiles.h5")
    ),
    y_pssm_h5_file: str = str(
        path.join(path.dirname(__file__), "../phosx/data/Y_PSSMs.h5")
    ),
    y_pssm_score_quantiles_h5_file: str = str(
        path.join(path.dirname(__file__), "../phosx/data/Y_PSSM_score_quantiles.h5")
    ),
    n_perm: int = 10000,
    s_t_n_top_kinases: int = 5,
    y_n_top_kinases: int = 10,
    min_n_hits: int = 4,
    s_t_min_quantile: float = 0.95,
    y_min_quantile: float = 0.90,
    no_upstream_activation_evidence: bool = False,
    metadata_h5_file: str = str(
        path.join(path.dirname(__file__), "../phosx/data/kinase_metadata.h5")
    ),
    a_loop_s_t_quantile_threshold: int = 0.95,
    a_loop_y_quantile_threshold: int = 0.90,
    upreg_redundancy_threshold: float = 0.5,
    downreg_redundancy_threshold: float = 0.5,
    decay_factor: float = 64,
    n_proc: int = 1,
    plot_figures: bool = False,
    out_plot_dir: str = "phosx_output",
    network_path=None,
    out_path=None,
):
    print("> Computing differential activity of Ser/Thr kinases...", file=sys.stderr)

    s_t_kinase_activity_df = pd.DataFrame()
    s_t_kinase_activity_df, s_t_assigned_substrates_df = compute_kinase_activities(
        seqrnk_file,
        s_t_pssm_h5_file,
        s_t_pssm_score_quantiles_h5_file,
        n_perm,
        s_t_n_top_kinases,
        min_n_hits,
        s_t_min_quantile,
        n_proc,
        plot_figures,
        out_plot_dir,
        True,
        False,
    )

    print("> Computing differential activity of Tyr kinases...", file=sys.stderr)

    y_kinase_activity_df = pd.DataFrame()
    y_kinase_activity_df, y_assigned_substrates_df = compute_kinase_activities(
        seqrnk_file,
        y_pssm_h5_file,
        y_pssm_score_quantiles_h5_file,
        n_perm,
        y_n_top_kinases,
        min_n_hits,
        y_min_quantile,
        n_proc,
        plot_figures,
        out_plot_dir,
        False,
        True,
    )

    activity_df = pd.concat([s_t_kinase_activity_df, y_kinase_activity_df], axis=0)

    if no_upstream_activation_evidence == False:
        print(
            "> Computing upstream activation evidence (upregulation) ...",
            file=sys.stderr,
        )

        upreg_activation_series, network = compute_activation_evidence(
            activity_df,
            s_t_assigned_substrates_df,
            y_assigned_substrates_df,
            s_t_pssm_h5_file,
            s_t_pssm_score_quantiles_h5_file,
            y_pssm_h5_file,
            y_pssm_score_quantiles_h5_file,
            metadata_h5_file,
            a_loop_s_t_quantile_threshold,
            a_loop_y_quantile_threshold,
            upreg_redundancy_threshold,
            decay_factor,
            True,
            n_proc,
            plot_figures,
            out_plot_dir,
            out_path,
        )

        print(
            "> Computing upstream activation evidence (downregulation) ...",
            file=sys.stderr,
        )

        downreg_activation_series, network = compute_activation_evidence(
            activity_df,
            s_t_assigned_substrates_df,
            y_assigned_substrates_df,
            s_t_pssm_h5_file,
            s_t_pssm_score_quantiles_h5_file,
            y_pssm_h5_file,
            y_pssm_score_quantiles_h5_file,
            metadata_h5_file,
            a_loop_s_t_quantile_threshold,
            a_loop_y_quantile_threshold,
            downreg_redundancy_threshold,
            decay_factor,
            False,
            n_proc,
            plot_figures,
            out_plot_dir,
            out_path,
        )

        activity_df = activity_df.rename(
            columns={"Activity Score": "Legacy Activity Score"}
        )
        activity_df.loc[upreg_activation_series.index, "Activity Score"] = (
            upreg_activation_series
        )
        activity_df.loc[downreg_activation_series.index, "Activity Score"] = (
            downreg_activation_series
        )

    # add the activity scores to the network
    network.add_activity_scores(activity_df["Activity Score"].to_dict())

    # export results
    if out_path == None:
        print(activity_df.to_csv(sep="\t", na_rep="NA", header=True, index=True))
    else:
        activity_df.to_csv(out_path, na_rep="NA", sep="\t", header=True, index=True)

    if network_path is not None:
        network_df = pd.DataFrame.from_dict(network.get_edges_dict(), orient="columns")
        network_df["complementarity"] = network_df["complementarity"].round(
            decimals=2
        )
        network_df.to_csv(
            f"{network_path}",
            sep="\t",
            na_rep="NA",
            header=True,
            index=False,
        )

    return activity_df


def main():
    print(
        f"""    
██████╗░██╗░░██╗░█████╗░░██████╗██╗░░██╗
██╔══██╗██║░░██║██╔══██╗██╔════╝╚██╗██╔╝
██████╔╝███████║██║░░██║╚█████╗░░╚███╔╝░
██╔═══╝░██╔══██║██║░░██║░╚═══██╗░██╔██╗░
██║░░░░░██║░░██║╚█████╔╝██████╔╝██╔╝╚██╗
╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═════╝░╚═╝░░╚═╝

Version 0.16.0
Copyright (C) 2025 Alessandro Lussana
Licence Apache 2.0

Command: {' '.join(sys.argv)}
    """,
        file=sys.stderr,
    )

    args = parse_phosx_args()

    phosx(
        args.seqrnk,
        args.s_t_pssm,
        args.s_t_pssm_quantiles,
        args.y_pssm,
        args.y_pssm_quantiles,
        args.n_permutations,
        args.s_t_n_top_kinases,
        args.y_n_top_kinases,
        args.min_n_hits,
        args.s_t_min_quantile,
        args.y_min_quantile,
        args.no_upstream_activation_evidence,
        args.kinase_metadata,
        args.a_loop_s_t_quantile_threshold,
        args.a_loop_y_quantile_threshold,
        args.upreg_redundancy_threshold,
        args.downreg_redundancy_threshold,
        args.decay_factor,
        args.n_proc,
        args.plot_figures,
        args.output_dir,
        args.network_path,
        args.output_path,
    )


if __name__ == "__main__":
    main()
