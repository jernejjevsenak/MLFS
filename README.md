<!-- README.md is generated from README.Rmd. Please edit that file -->

# MLFS

*Machine Learning Forest Simulator (MLFS)* is a data-driven forest development model implemented as an R package. MLFS aims to remove the need for hand-tuned model parametrization and provide an easy-to-use, freely available tool applicable across forest types—from even-aged monocultures to mixed, vertically structured stands.

MLFS is age- and spatially independent. It uses repeated forest inventory data (at least two periods) to:
1) model basal area increment (BAI),
2) update tree and crown heights at each simulation step,
3) simulate individual-tree mortality,
4) simulate ingrowth of new trees, and
5) represent harvesting.

The main inputs are forest inventory plot data and site descriptors, which are used to train modular sub-models.

## Key features

- **Data-driven**: learns directly from repeated inventory data.
- **Modular workflow**: separate sub-models for growth, height dynamics, mortality, ingrowth, and harvesting.
- **Broad applicability**: works from even-aged to mixed, multi-layered forests.
- **R-native**: integrates with common R tooling for analysis and reproducibility.

## Installation

Install the released version from CRAN:

```r
install.packages("MLFS")
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("jernejjevsenak/MLFS")
```

## Input data requirements (at a glance)

- Forest inventory plots with at least two measurement periods.
- Tree-level attributes (e.g., species, DBH; optional heights/crown metrics if available).
- Site descriptors used for training sub-models.

## Getting started

After installation:

```r
library(MLFS)

# Explore package documentation
help(package = "MLFS")
```

## A typical workflow is:

- Prepare and validate inventory and site tables.
- Train sub-models using repeated measurements.
- Run simulations over chosen time steps.
- Summarize and visualize outputs at tree/plot/stand levels.

## Help & issues

- Report bugs or request features via GitHub Issues: https://github.com/jernejjevsenak/MLFS/issues

## How to cite

```r
citation("MLFS")
```

## Authors

Jernej Jevšenak

