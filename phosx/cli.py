#!/usr/bin/env python3

# Alessandro Lussana <alussana@ebi.ac.uk>

import phosx
import argparse
from os import path

def parse_phosx_args():
    parser = argparse.ArgumentParser(
                prog='phosx',
                description='Kinase activity inference from phosphosproteomics data based on substrate sequence specificity',
                epilog=''
    )    
    parser.add_argument(
        'seqrnk',
        type=str,
        help='Path to the seqrnk file.'
    )
    parser.add_argument(
        '-p', '--pssm',
        type=str,
        default=str(path.join(path.dirname(__file__), 'data/pssm.h5')),
        help='Path to the h5 file storing custom PSSMs; defaults to built-in PSSMs'
    )
    parser.add_argument(
        '-n', '--n_permutations',
        type=int,
        default=10000,
        help='Number of random permutations; defaults to 10000'
    )
    parser.add_argument(
        '-c', '--n_proc',
        type=int,
        default=1,
        help='Number of cores used for multithreading; defaults to 1'
    )
    parser.add_argument(
        '--plot_figures',
        action='store_true',
        help='Output files prefix; defaults to \'phosx_output/\''
    )
    parser.add_argument(
        '-o', '--output_prefix',
        type=str,
        default='phosx_output/',
        help='Output files prefix; only relevant if used with --plot_figures; defaults to \'phosx_output/\''
    )
    args = parser.parse_args()
    return args

def main():
    args = parse_phosx_args()
    phosx.compute_kinase_activity(
        args.seqrnk,
        args.pssm,
        args.n_permutations,
        args.n_proc,
        args.plot_figures,
        args.output_prefix
    )

if __name__ == '__main__':
    main()