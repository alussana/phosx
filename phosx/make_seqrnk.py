#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import requests
import warnings
from functools import lru_cache


@lru_cache
def get_fasta(uniprot_ac: str):
    url_prefix = "https://rest.uniprot.org/uniprotkb/"
    r = requests.get(f"{url_prefix}{uniprot_ac}.fasta")
    fasta = r.content.decode()
    lines = fasta.split("\n")
    lines.pop(0)
    sequence = "".join(lines)
    return sequence


def get_phosphosite(uniprot_ac: str, pos: int, before: int, after: int):
    sequence = get_fasta(uniprot_ac)
    seq_length = len(sequence)
    before_padding = "_" * before
    after_padding = "_" * after
    sequence = before_padding + sequence + after_padding
    try:
        phosphosite = sequence[pos - 1 : pos + before + after]
        return phosphosite
    except:
        warnings.warn(
            f"""
            Phosphosite position exceedes sequence length.
            pos: {pos}
            sequence length: {seq_length}
            """
        )


def make_seqrnk_entry(series: pd.Series) -> pd.Series:
    uniprot_ac = series["UniProtAC"]
    pos = series["Position"]
    sequence = get_phosphosite(uniprot_ac=uniprot_ac, pos=pos, before=5, after=4)
    score = series["Score"]
    seqrnk_entry_series = pd.Series(
        {"Sequence": sequence, "Score": score}, name=series.name
    )
    return seqrnk_entry_series


def parse_args():
    parser = argparse.ArgumentParser(
        prog="make-seqrnk",
        description="Make a seqrnk file to be used to compute differential kinase activity with PhosX",
        epilog="",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=sys.stdin,
        help="Path of the input phosphosites to be converted in seqrnk format. It should be a TSV file where the 1st column is the UniProtAC (str), the 2nd is the sequence coordinate (int), and the 3rd is the logFC (float); defaults to STDIN",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=sys.stdout,
        help="Path of the seqrnk file; if not specified it will be printed in STDOUT",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    df = pd.read_csv(filepath_or_buffer=args.input, sep="\t", header=None)
    df.columns = ["UniProtAC", "Position", "Score"]
    df["Position"] = df["Position"].astype(int)
    seqrnk = df.apply(make_seqrnk_entry, axis=1)
    seqrnk = seqrnk.sort_values("Score", ascending=False)
    seqrnk.to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
