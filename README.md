# PhosX

Kinase activity inference from phosphosproteomics data based on substrate sequence specificity.

## Installation

### Build from source (requires [Poetry](https://python-poetry.org))

```
poetry build
pip install dist/*.whl
```

### Install from [PyPI](https://pypi.org)

```
pip install phosx
```

## Usage

```
phosx [-h] [-v] [-p PSSM] [-n N_PERMUTATIONS] [-c N_PROC] [--plot-figures] [-o OUTPUT_DIR] seqrnk
```
```
positional arguments:
  seqrnk                Path to the seqrnk file.

options:
  -h, --help            show this help message and exit
  -v, --version         Print package version and exit
  -p PSSM, --pssm PSSM  Path to the h5 file storing custom PSSMs; defaults to built-in PSSMs
  -n N_PERMUTATIONS, --n-permutations N_PERMUTATIONS
                        Number of random permutations; defaults to 10000
  -c N_PROC, --n-proc N_PROC
                        Number of cores used for multithreading; defaults to 1
  --plot-figures        Save figures in svg format; see also --output_dir
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output files directory; only relevant if used with --plot_figures; defaults to 'phosx_output/'
```

## TODO

- [ ] **Filter temporary warnings** - when plotting figures from `phosx.pssm_enrichments.ks_statistic` and `phosx.pssm_enrichments.compute_ks_pvalues`, many warnings will be thrown (`UserWarning: The figure layout has changed to tight`). This is a Matplotlib known issue: [https://github.com/matplotlib/matplotlib/issues/26290](https://github.com/matplotlib/matplotlib/issues/26290)
- [x] **Create dirs** - if not existent, for output files
- [ ] **Add docs**
- [ ] **FDR** - Report also the FDR q-value in the output table  