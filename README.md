<p align="center">
  <img width="190" src="https://i.imgur.com/OzGTvkt.png">
  <br>
  Kinase activity inference from phosphosproteomics data based on substrate sequence specificity
  <br><br>
</p>

![Build and publish to PyPI badge](https://github.com/alussana/phosx/actions/workflows/build-and-publish-to-pypi.yml/badge.svg)

> Current version: `0.7.0`

> Research paper: [https://doi.org/10.1101/2024.03.22.586304](https://doi.org/10.1101/2024.03.22.586304)

> Benchmark repository: [https://github.com/alussana/phosx-benchmark](https://github.com/alussana/phosx-benchmark)

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

```bash
phosx [-h] [-yp Y_PSSM] [-stp S_T_PSSM] [-yq Y_PSSM_QUANTILES] [-stq ST_PSSM_QUANTILES] [-n N_PERMUTATIONS] [-stk S_T_N_TOP_KINASES] [-yk Y_N_TOP_KINASES] [-m MIN_N_HITS] [-c N_PROC] [--plot-figures] [-d OUTPUT_DIR] [-o OUTPUT_PATH] [-v] seqrnk
```
```bash
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
  -stq ST_PSSM_QUANTILES, --st-pssm-quantiles ST_PSSM_QUANTILES
                        Path to the h5 file storing custom Ser/Thr kinases PSSM score quantile distributions under the key 'pssm_scores'; defaults to built-in PSSM scores quantiles
  -n N_PERMUTATIONS, --n-permutations N_PERMUTATIONS
                        Number of random permutations; default: 1000
  -stk S_T_N_TOP_KINASES, --s-t-n-top-kinases S_T_N_TOP_KINASES
                        Number of top-scoring Ser/Thr kinases potentially associatiated to a given phosphosite; default: 5
  -yk Y_N_TOP_KINASES, --y-n-top-kinases Y_N_TOP_KINASES
                        Number of top-scoring Tyr kinases potentially associatiated to a given phosphosite; default: 5
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

# Cite

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