#' An example of joined national forest inventory and site data that is used
#' within the MLFS
#'
#' This is simulated data that reassemble the national forest inventory data.
#' We use it to show how to run examples for simulating harvesting.
#'
#' @format A data frame with 5949 rows and 10 variables:
#' \describe{
#'   \item{species}{species name}
#'   \item{year}{year in which plot was visited}
#'   \item{plotID}{a unique identifier for plot}
#'   \item{treeID}{a unique identifier for tree}
#'   \item{speciesGroup}{identifier for species group}
#'   \item{code}{status of a tree: 0 (normal), 1(harvested), 2(dead), 3 (ingrowth)}
#'   \item{volume_mid}{tree volume in m3 in the middle of a simulation step}
#'   \item{weight_mid}{upscale weight to calculate hectare values in the middle of a simulation step}
#'   \item{BA_mid}{basal area of individual trees in m2 in the middle of a simulation step}
#'   \item{protected}{logical, 1 if protected, otherwise 0}
#' }
#' @export
"data_v5"
