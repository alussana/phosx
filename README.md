<p align="center">
  <img width="230" src="https://i.imgur.com/OzGTvkt.png">
  <br>
  Kinase activity inference from phosphosproteomics data based on substrate sequence specificity
  <br><br>
</p>

![Build and publish to PyPI badge](https://github.com/alussana/phosx/actions/workflows/build-and-publish-to-pypi.yml/badge.svg)

> Current version: `v0.4.1`
> 
> NOTE: this software is still in development. Unexpected behaviors may occur.


# Installation

## From [PyPI](https://pypi.org)

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
phosx [-h] [-p PSSM] [-q PSSM_QUANTILES] [-n N_PERMUTATIONS] [-k N_TOP_KINASES] [-c N_PROC] [--plot-figures] [-o OUTPUT_DIR] [-v] seqrnk
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
                        Number of random permutations; defaults to 10000
  -k N_TOP_KINASES, --n-top-kinases N_TOP_KINASES
                        Number of top-scoring kinases potentially associatiated to a given phosphosite; defaults to 15
  -m MIN_N_HITS, --min-n-hits MIN_N_HITS
                        Minimum number of phosphosites associated with a kinase for the kinase to be considered in the analysis; defaults to 15
  -c N_PROC, --n-proc N_PROC
                        Number of cores used for multithreading; defaults to 1
  --plot-figures        Save figures in svg format; see also --output_dir
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output files directory; only relevant if used with --plot_figures; defaults to 'phosx_output/'
  -v, --version         Print package version and exit
```

# TODO

- [ ] **Filter temporary warnings** - when plotting figures from `phosx.pssm_enrichments.ks_statistic` and `phosx.pssm_enrichments.compute_ks_pvalues`, many warnings will be thrown (`UserWarning: The figure layout has changed to tight`). This is a Matplotlib known issue: [https://github.com/matplotlib/matplotlib/issues/26290](https://github.com/matplotlib/matplotlib/issues/26290)
- [ ] **Add docs**
- [ ] **tqdm** - Add optional progress bar for permutations
- [x] **Null KS** - Report also null KS mean in the output table
- [ ] **FDR** - Report also the FDR q-value in the output table