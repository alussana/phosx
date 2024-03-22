<p align="center">
  <img width="190" src="https://i.imgur.com/OzGTvkt.png">
  <br>
  Kinase activity inference from phosphosproteomics data based on substrate sequence specificity
  <br><br>
</p>

![Build and publish to PyPI badge](https://github.com/alussana/phosx/actions/workflows/build-and-publish-to-pypi.yml/badge.svg)

> Current version: `v0.5.3`

> Research paper: [ ]

# Overview

<p align="center">
<br>
  <img width="900" src="https://i.imgur.com/6DdMDom.png">
  <br>
</p>

PhosX infers differential kinase activities from phosphoproteomics data without requiring any prior knowledge database of kinase-phosphosite associations. PhosX assigns the detected phosphopeptides to potential upstream kinases based on experimentally determined substrate sequence specificities, and it tests the enrichment of a kinase's potential substrates in the extremes of a ranked list of phosphopeptides using a Kolmogorov-Smirnov-like statistic. A _p_ value for this statistic is extracted empirically by random permutations of the phosphosite ranks.

# Installation

## From [PyPI](https://pypi.org/project/phosx/)

```
pip install phosx
```

## From source (requires [Poetry](https://python-poetry.org))

```
poetry build
pip install dist/*.whl
```

# Usage

```
phosx [-h] [-p PSSM] [-q PSSM_QUANTILES] [-n N_PERMUTATIONS] [-k N_TOP_KINASES] [-m MIN_N_HITS] [-c N_PROC] [--plot-figures] [-d OUTPUT_DIR] [-o OUTPUT_PATH] [-v] seqrnk
```
```
positional arguments:
  seqrnk                Path to the seqrnk file.

options:
  -h, --help            show this help message and exit
  -p PSSM, --pssm PSSM  Path to the h5 file storing custom PSSMs; defaults to built-in PSSMs
  -q PSSM_QUANTILES, --pssm-quantiles PSSM_QUANTILES
                        Path to the h5 file storing custom PSSM score quantile distributions under the key 'pssm_scores'; defaults to built-in PSSM scores quantiles
  -n N_PERMUTATIONS, --n-permutations N_PERMUTATIONS
                        Number of random permutations; default: 1000
  -k N_TOP_KINASES, --n-top-kinases N_TOP_KINASES
                        Number of top-scoring kinases potentially associatiated to a given phosphosite; default: 8
  -m MIN_N_HITS, --min-n-hits MIN_N_HITS
                        Minimum number of phosphosites associated with a kinase for the kinase to be considered in the analysis; default: 4
  -c N_PROC, --n-proc N_PROC
                        Number of cores used for multithreading; default: 1
  --plot-figures        Save figures in pdf format; see also --output_dir
  -d OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output files directory; only relevant if used with --plot_figures; defaults to 'phosx_output/'
  -o OUTPUT_PATH, --output-path OUTPUT_PATH
                        Main output table; if not specified it will be printed in STDOUT
  -v, --version         Print package version and exit
```

Minimal example to run PhosX with default parameters on an example dataset, using up to 8 cores, and redirecting the output table to `kinase_activities.tsv`:

```bash
phosx -c 8 tests/seqrnk/koksal2018_log2.fold.change.8min.seqrnk > kinase_activities.tsv
```