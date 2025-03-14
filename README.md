<p align="center">
  <img width="190" src="https://i.imgur.com/OzGTvkt.png">
  <br>
  Kinase activity inference from phosphosproteomics data based on substrate sequence specificity
  <br><br>
</p>

![Build and publish to PyPI badge](https://github.com/alussana/phosx/actions/workflows/build-and-publish-to-pypi.yml/badge.svg)

> Current version: `0.9.2`

> Research paper: [https://doi.org/10.1101/2024.03.22.586304](https://doi.org/10.1101/2024.03.22.586304)

> Benchmark: [https://github.com/alussana/phosx-benchmark](https://github.com/alussana/phosx-benchmark)

> Data: [https://github.com/alussana/kinase_pssms](https://github.com/alussana/kinase_pssms)

# Overview

<p align="center">
<br>
  <img width="900" src="https://i.imgur.com/6DdMDom.png">
  <br>
</p>

PhosX infers differential kinase activities from phosphoproteomics data without requiring any prior knowledge database of kinase-phosphosite associations. PhosX assigns the detected phosphopeptides to potential upstream kinases based on experimentally determined substrate sequence specificities, and it tests the enrichment of a kinase's potential substrates in the extremes of a ranked list of phosphopeptides using a Kolmogorov-Smirnov-like statistic. A _p_ value for this statistic is extracted empirically by random permutations of the phosphosite ranks.

# Installation

## From [PyPI](https://pypi.org/project/phosx/)

```bash
pip install phosx
```

## From source (requires [Poetry](https://python-poetry.org))

```
poetry build
pip install dist/*.whl
```

# Usage

PhosX can be used as a command line tool (`phosx`) with minimal effort. Its output is redirected by default in `STDOUT`, making it easy to use in bioinformatics pipelines. Alternatively, the user can specify an output filename (option `-o`). 

Example: run PhosX with default parameters on an example dataset, using up to 8 cores, and redirecting the output table to `kinase_activities.tsv`:

```bash
phosx -c 8 tests/seqrnk/koksal2018_log2.fold.change.8min.seqrnk > kinase_activities.tmp
```

A full description of the command line options can be viewed with `phosx -h`:

```bash 
usage: phosx [-h] [-yp Y_PSSM] [-stp S_T_PSSM] [-yq Y_PSSM_QUANTILES] [-stq S_T_PSSM_QUANTILES] [-n N_PERMUTATIONS] [-stk S_T_N_TOP_KINASES] [-yk Y_N_TOP_KINASES] [-mh MIN_N_HITS] [-mp MIN_QUANTILE] [-c N_PROC] [--plot-figures] [-d OUTPUT_DIR]
             [-o OUTPUT_PATH] [-v]
             seqrnk

Kinase activity inference from phosphosproteomics data based on substrate sequence specificity

positional arguments:
  seqrnk                Path to the seqrnk file.

options:
  -h, --help            show this help message and exit
  -yp Y_PSSM, --y-pssm Y_PSSM
                        Path to the h5 file storing custom Tyr PSSMs; defaults to built-in PSSMs
  -stp S_T_PSSM, --s-t-pssm S_T_PSSM
                        Path to the h5 file storing custom Ser/Thr PSSMs; defaults to built-in PSSMs
  -yq Y_PSSM_QUANTILES, --y-pssm-quantiles Y_PSSM_QUANTILES
                        Path to the h5 file storing custom Tyr kinases PSSM score quantile distributions under the key 'pssm_scores'; defaults to built-in PSSM scores quantiles
  -stq S_T_PSSM_QUANTILES, --s-t-pssm-quantiles S_T_PSSM_QUANTILES
                        Path to the h5 file storing custom Ser/Thr kinases PSSM score quantile distributions under the key 'pssm_scores'; defaults to built-in PSSM scores quantiles
  -n N_PERMUTATIONS, --n-permutations N_PERMUTATIONS
                        Number of random permutations; default: 1000
  -stk S_T_N_TOP_KINASES, --s-t-n-top-kinases S_T_N_TOP_KINASES
                        Number of top-scoring Ser/Thr kinases potentially associatiated to a given phosphosite; default: 5
  -yk Y_N_TOP_KINASES, --y-n-top-kinases Y_N_TOP_KINASES
                        Number of top-scoring Tyr kinases potentially associatiated to a given phosphosite; default: 5
  -mh MIN_N_HITS, --min-n-hits MIN_N_HITS
                        Minimum number of phosphosites associated with a kinase for the kinase to be considered in the analysis; default: 4
  -mp MIN_QUANTILE, --min-quantile MIN_QUANTILE
                        Minimum PSSM score quantile that a phosphosite has to satisfy to be potentially assigned to a kinase; default: 0.95
  -c N_PROC, --n-proc N_PROC
                        Number of cores used for multithreading; default: 1
  --plot-figures        Save figures in pdf format; see also --output_dir
  -d OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output files directory; only relevant if used with --plot_figures; defaults to 'phosx_output/'
  -o OUTPUT_PATH, --output-path OUTPUT_PATH
                        Main output table; if not specified it will be printed in STDOUT
  -v, --version         Print package version and exit
```

For a full description of the method please see the [Method](#method) section and the [manuscript](https://doi.org/10.1101/2024.03.22.586304).

## Input

### _seqrnk_

PhosX's input format is a simple text file that we name _seqrnk_. It consists of 2 tab-separated columns containing phosphopeptide sequences and values, respectively. The values should be biologically relevant measures of differential phosphorylation, typically intensity log fold changes as obtained when comparing two conditions in mass spectrometry experiments. Amino acid sequences should be of length $10$, with the phosphorylated residue in position $6$ (1-based), in order to match the Position Specific Scoring Matrix (PSSM) models (see [Ser/Thr PSSMs](https://doi.org/10.1038/s41586-022-05575-3) & [Tyr PSSMs](https://doi.org/10.1038/s41586-024-07407-y)). Undefined amino acids are represented by the character `_`. Every other residue is represented by the corresponding 1-letter symbol according to the IUPAC nomenclature for amino acids and additional phosphorylated Serine, Threonine or Tyrosine residues are represented with the symbols `s`, `t`, and `y`, respectively. Phosphorylated residues that act as potential priming sites and are therefore not in the $6^{th}$ position of the peptide are represented with lowercase letters. An example is included in this repository:

```bash
$ head tests/seqrnk/koksal2018_log2.fold.change.8min.seqrnk

QEEAEYVRAL      5.646644
ANFSAYPSEE      4.33437
YLNRNYWEKK      4.174151
AENAEYLRVA      3.685413
STYTSYPKAE      3.491975
SFLQRYSSDP      3.295341
AAEPGSPTAA      3.202242
EPAHAYAQPQ      3.160899
RQKSTYTSYP      3.114077
ETKSLYPSSE      3.04653
```

Alongside the main program, this package also installs `make-seqrnk`. This utily can be used to help generating a _seqrnk_ file given a list of phosphosites, each one identified by a UniProt Acession Number and residue coordinate. `make-seqrnk` will query the [UniProt](https://www.uniprot.org) database to fetch the appropriate subsequences to build the _seqrnk_ file. See `make-seqrnk -h` for more details: 

```bash
usage: make-seqrnk [-h] [-i INPUT] [-o OUTPUT]

Make a seqrnk file to be used to compute differential kinase activity with PhosX

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path of the input phosphosites to be converted in seqrnk format. It should be a TSV file where the 1st column is the UniProtAC (str), the 2nd is the sequence coordinate (int), and the 3rd is the logFC (float); defaults to STDIN
  -o OUTPUT, --output OUTPUT
                        Path of the seqrnk file; if not specified it will be printed in STDOUT
```

Run an example:

```bash
cat tests/p_list/15_3.tsv | make-seqrnk > 15_3.seqrnk
```

### PSSMs

PhosX estimates the affinity between human kinases and phosphopeptides based on the substrate sequence specificity encoded in Position Specific Scoring Matrices (PSSMs). A kinase PSSM is a $10 \times 23$ matrix containing amino acid affinity scores at each one of the $10$ positions of the substrate. The 6\textsuperscript{th} position corresponds to the modified residue and should have non-$0$ values only for Serine, Threonine (for Ser/Thr kinases), or Tyrosine (for Tyr kinases) residues.

PhosX comes with built-in, default PSSMs for human kinases, that can be found at `phosx/data/*_PSSMs.h5`. The user can also run PhosX using custom PSSMs, whose path can be specified with the options `-yp` and `-stp`. 

For convenience, here is a function that can be used to open and inspect the structure of the Hierarchical Data Format version 5 (HDF5) files where the PSSMs are stored:

```python
import h5py


AA_LIST = [
  "G","P","A","V","L","I",
  "M","C","F","Y","W","H",
  "K","R","Q","N","E","D",
  "S","T","s","t","y",
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
```

Similarly, PhosX also has built-in kinase PSSM scores quantile distributions computed on a reference human phosphoproteome from the [PhosphositePlus](https://www.phosphosite.org/) database. These can be found at `phosx/data/*_PSSM_score_quantiles.h5`. When supplying custom PSSMs, it is necessary to also specify the appropriate background distributions with the options `-yq` and `-stq`. Inspect the HDF5 files containing the background PSSM scores with:

```python
import h5py


def read_pssm_score_quantiles(pssm_score_quantiles_h5_file: str):

    pssm_bg_scores_df = pd.read_hdf(pssm_score_quantiles_h5_file, key="pssm_scores")

    return pssm_bg_scores_df
```

## Output

PhosX's main output is a text file reporting the computed kinase activities with associated statistics as described in the [Method](#method) section. For each kinase, the KS	statistics, the _p_ value, the FDR _q_ value, and the	Activity Score are reported. Note that kinases for which an Activity Score could not be computed (for example for lack of matching phosphopeptides in the input data) are omitted in the current version. See an output example from the command executed [above](#usage):

```bash
$ head kinase_activities.tmp

        KS      p value FDR q value     Activity Score
AAK1    -0.2476131530554456     0.533   1.0     -0.2732727909734277
ACVR1B  -0.36580307230946174    0.078   1.0     -1.1079053973095196
ACVR2A  0.2259224236806207      0.439   1.0     0.35753547975787864
ACVR2B  0.47014516632215597     0.019   1.0     1.7212463990471711
ALK2    -0.25190558854944195    0.276   1.0     -0.5590909179347823
ALPHAK3 -0.2875279855211264     0.358   1.0     -0.44611697335612566
ALPK3   0.5759513630398431      0.041   1.0     1.3872161432802645
AMPKA2  0.4401107718873606      0.299   1.0     0.5243288116755703
ATM     -0.49337471491068263    0.096   1.0     -1.0177287669604316
```

Additionally, PhosX can also save plots of the of the weighted running sum and of the KS statistic compared to its empirical null distribution, similarly to the ones show [above](#overview), for each kinase. To enable this behavior the option `--plot-figures` must be specified. A custom directory to save the plots can be passed with `-d`.

# Method

## Phosphopeptide scoring

For each kinase PSSM, a score is assigned to each phosphopeptide sequence $S$ that quantifies its similarity to the PSSM. First, a "raw PSSM score" is computed as:

```math
\texttt{score}(S,k) := \prod_{i=-5}^{4}  
\begin{cases}
    M^{k}_{i,S_i}, & \text{if } S_i \neq \texttt{'\_'} \\
    1, & \text{if } S_i = \texttt{'\_'} 
\end{cases}
```

where $S_i$ is the amino acid residue at position $i$ of the phosphopeptide sequence $S$; $M^k_{i,j}$ is the value of the PSSM for kinase $k$ at position $i$ for residue $j$. Raw PSSM scores for each kinase are then transformed between $0$ and $1$ based on the quantile they fall in, considering a background distribution of proteome-wide raw PSSM scores. For each kinase, phosphopeptides with raw PSSM score equal to $0$ are discarded, and the remaining are used to determine the values of the $10,000$-quantiles of the raw PSSM score distribution. The background $10,000$-quantiles raw PSSM scores for each kinase PSSM are used to derive the final PSSM scores for each phosphopeptide.
￼
## Weighted running sum statistic

PhosX uses the PSSM scores to link kinases to their potential substrates. Each phosphopeptide is assigned as potential target to its $n$ top-scoring kinases, with default value of $10$. With the method has little sensitivity to this parameter in the range $[5,15]$. The activity change of a given kinase is estimated by calculating a running sum statistic over the ranked list of phosphosites, and by estimating its significance based on an empirical distribution generated by random permutations of the ranks. Let $C$ be the set of indexes corresponding to the ranked phosphosites associated with kinase $k$; $N$ the total number of phosphosites; $N_h$ the size of $C$; $r_i$ the value of the ranking metric of the phosphosite at rank $i$, where $r_0$ is the highest value. Then, the running sum ($RS$) up to the phosphosite rank $n$ is given by

```math
RS(k,n) := \sum_{i=0}^{n}  
\begin{cases}
    \frac{|r_i|}{N_R}, & \text{if } i \in C \\
    -\frac{1}{N - N_h}, & \text{if } i \not\in C 
\end{cases}
```

where

```math
N_R = \sum_{i \in C} |r_i|
```

The kinase enrichment score ($ES$) corresponds to the maximum deviation from $0$ of $RS$.

## Permutation tests

For each kinase, PhosX computes an empirical _p_ value of the $ES$ by generating a null distribution of the $ES$ through random permutations of the phosphosite ranks. A False Discovery Rate (FDR) _q_ value is also calculated by applying the Bonferroni method considering the number of kinases independently tested. The number of permutations is a tunable parameter but we recommend performing at least $10^4$ random permutations to be able to compute FDR values $< 0.05$.

## Kinase activity score

The activity score for a given kinase is defined as:

```math
Activity = -\log_{10}{(p)} * \texttt{sign}(ES)
```

where \texttt{sign} is the sign function, and $-log_{10}{(p)}$ is capped at the smallest computable \textit{p} value different from $0$, _i.e._ the inverse of the number of random permutations.
Activity scores greater than $0$ denote kinase activation, while the opposite corresponds to kinase inhibition.

# Cite

Please cite one of the following references if you use PhosX in your work.

## BioRxiv

BibTeX:

```bibtex
@article{Lussana2024,
  title = {PhosX: data-driven kinase activity inference from phosphoproteomics experiments},
  url = {http://dx.doi.org/10.1101/2024.03.22.586304},
  DOI = {10.1101/2024.03.22.586304},
  publisher = {Cold Spring Harbor Laboratory},
  author = {Lussana,  Alessandro and Petsalaki,  Evangelia},
  year = {2024},
  month = mar 
}
```

## Bioinformatics

BibTeX:

```bibtex
@article{10.1093/bioinformatics/btae697,
    author = {Lussana, Alessandro and Müller-Dott, Sophia and Saez-Rodriguez, Julio and Petsalaki, Evangelia},
    title = {PhosX: data-driven kinase activity inference from phosphoproteomics experiments},
    journal = {Bioinformatics},
    volume = {40},
    number = {12},
    pages = {btae697},
    year = {2024},
    month = {11},
    abstract = {The inference of kinase activity from phosphoproteomics data can point to causal mechanisms driving signalling processes and potential drug targets. Identifying the kinases whose change in activity explains the observed phosphorylation profiles, however, remains challenging, and constrained by the manually curated knowledge of kinase–substrate associations. Recently, experimentally determined substrate sequence specificities of human kinases have become available, but robust methods to exploit this new data for kinase activity inference are still missing. We present PhosX, a method to estimate differential kinase activity from phosphoproteomics data that combines state-of-the-art statistics in enrichment analysis with kinases’ substrate sequence specificity information. Using a large phosphoproteomics dataset with known differentially regulated kinases we show that our method identifies upregulated and downregulated kinases by only relying on the input phosphopeptides’ sequences and intensity changes. We find that PhosX outperforms the currently available approach for the same task, and performs better or similarly to state-of-the-art methods that rely on previously known kinase–substrate associations. We therefore recommend its use for data-driven kinase activity inference.PhosX is implemented in Python, open-source under the Apache-2.0 licence, and distributed on the Python Package Index. The code is available on GitHub (https://github.com/alussana/phosx).},
    issn = {1367-4811},
    doi = {10.1093/bioinformatics/btae697},
    url = {https://doi.org/10.1093/bioinformatics/btae697},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/40/12/btae697/60972735/btae697.pdf},
}
```