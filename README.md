<p align="center">
  <img width="190" src="https://i.imgur.com/OzGTvkt.png">
  <br>
  Kinase activity inference from phosphosproteomics data based on substrate sequence specificity
  <br><br>
</p>

![Build and publish to PyPI badge](https://github.com/alussana/phosx/actions/workflows/build-and-publish-to-pypi.yml/badge.svg)

> Current version: `0.8.0`

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

Run PhosX with default parameters on an example dataset, using up to 8 cores, and redirecting the output table to `kinase_activities.tsv`:

```bash
phosx -c 8 tests/seqrnk/koksal2018_log2.fold.change.8min.seqrnk > kinase_activities.tsv
```

See the full list of command line options with `phosx -h`.

Alongside the main program, this package also installs `make-seqrnk`. This utily can be used to easily generate a _seqrnk_ file, which is used as input by PhosX, given a list of phosphosites, each one identified by a UniProtAC and residue coordinate. `make-seqrnk` will query the UniProt database to fetch the appropriate subsequences and build the _seqrnk_ file. Run `make-seqrnk -h` for more details, or see an example with: 

```bash
cat tests/p_list/15_3.tsv | make-seqrnk > 15_3.seqrnk
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