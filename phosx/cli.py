#!/usr/bin/env python3

import argparse
from os import path
import sys
from phosx.pssm_enrichments import kinase_activities


def parse_phosx_args():
    parser = argparse.ArgumentParser(
        prog="phosx",
        description="Kinase activity inference from phosphosproteomics data based on substrate sequence specificity",
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
        "-n",
        "--n-permutations",
        type=int,
        default=1000,
        help="Number of random permutations; default: 1000",
    )
    parser.add_argument(
        "-stk",
        "--s-t-n-top-kinases",
        type=int,
        default=5,
        help="Number of top-scoring Ser/Thr kinases potentially associatiated to a given phosphosite; default: 5",
    )
    parser.add_argument(
        "-yk",
        "--y-n-top-kinases",
        type=int,
        default=5,
        help="Number of top-scoring Tyr kinases potentially associatiated to a given phosphosite; default: 5",
    )
    parser.add_argument(
        "-m",
        "--min-n-hits",
        type=int,
        default=4,
        help="Minimum number of phosphosites associated with a kinase for the kinase to be considered in the analysis; default: 4",
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
        "-o",
        "--output-path",
        type=str,
        default=None,
        help="Main output table; if not specified it will be printed in STDOUT",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="v0.7.0",
        help="Print package version and exit",
    )
    args = parser.parse_args()
    return args


def main():
    print(
        f"""    
  ██████╗░██╗░░██╗░█████╗░░██████╗██╗░░██╗
  ██╔══██╗██║░░██║██╔══██╗██╔════╝╚██╗██╔╝
  ██████╔╝███████║██║░░██║╚█████╗░░╚███╔╝░
  ██╔═══╝░██╔══██║██║░░██║░╚═══██╗░██╔██╗░
  ██║░░░░░██║░░██║╚█████╔╝██████╔╝██╔╝╚██╗
  ╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═════╝░╚═╝░░╚═╝
  
  Version 0.7.0
  Copyright (C) 2024 Alessandro Lussana
  Licence Apache 2.0
  
  Command: {' '.join(sys.argv)}
    """,
        file=sys.stderr,
    )

    args = parse_phosx_args()
    kinase_activities(
        args.seqrnk,
        args.s_t_pssm,
        args.s_t_pssm_quantiles,
        args.y_pssm,
        args.y_pssm_quantiles,
        args.n_permutations,
        args.s_t_n_top_kinases,
        args.y_n_top_kinases,
        args.min_n_hits,
        args.n_proc,
        args.plot_figures,
        args.output_dir,
        args.output_path,
    )


if __name__ == "__main__":
    main()
