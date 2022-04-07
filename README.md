
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MLFS

Machine Learning Forest Simulator (MLFS) is the first complete
data-driven forest development model, organized as R package. The main
motivation behind the development was to remove the need for model
parametrization and to provide easy to use and freely available tool,
applicable to all forest ecosystems, from even-aged monocultures to
mixed forests with diverse vertical structures. MLFS is freely
available, age- and spatially-independent forest development tool. It
requires data from at least two inventory periods, which is used to 1)
model basal area increments (BAI), 2) update tree and crown heights in
each simulation step, 3) simulate individual tree mortality, 4) ingrowth
of new trees, and 5) harvesting. The main input tables consist of forest
inventory plot data, and site descriptors, which are used to train
specific sub-models.

## Installation

You can install MLFS using:

``` r
library("devtools")
devtools::install_github("jernejjevsenak/MLFS") # current version under development

install.packages("MLFS") # from CRAN
```

## Authors

-   **Jernej Jevšenak**
