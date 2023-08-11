# PhosX

Kinase activity inference from phosphosproteomics data based on substrate sequence specificity

## Installation

### Build from source (requires [Poetry](https://python-poetry.org))

```
poetry build
pip install dist/*.whl
```

### Install from PyPI (currently live only on [Test PyPI](https://test.pypi.org))

```
pip install phosx
```

## Usage

```
phosx [-h] [-v] [-p PSSM] [-n N_PERMUTATIONS] [-c N_PROC] [--plot_figures] [-o OUTPUT_DIR] seqrnk
```

## TODO

- [ ] **Filter temporary warnings** - when plotting figures from `phosx.pssm_enrichments.ks_statistic` and `phosx.pssm_enrichments.compute_ks_pvalues`, many warnings will be thrown (`UserWarning: The figure layout has changed to tight`). This is a Matplotlib known issue: [https://github.com/matplotlib/matplotlib/issues/26290](https://github.com/matplotlib/matplotlib/issues/26290)
- [x] **Create dirs** - if not existent, for output files
- [ ] **Add docs**
- [ ] **FDR** - Report also the FDR q-value in the output table  