#' An example of joined national forest inventory and site data that is used
#' within the MLFS
#'
#' This is simulated data that reassemble the national forest inventory data.
#' We use it to show how to run examples for tree and crown height predictions.
#' The difference between data_v2 and data_v3 is that in data_v3, tree heights
#' are already predicted
#'
#' @format A data frame with 6948 rows and 14 variables:
#' \describe{
#'   \item{plotID}{a unique identifier for plot}
#'   \item{treeID}{a unique identifier for tree}
#'   \item{year}{year in which plot was visited}
#'   \item{speciesGroup}{identifier for species group}
#'   \item{code}{status of a tree: 0 (normal), 1(harvested), 2(dead), 3 (ingrowth)}
#'   \item{species}{species name}
#'   \item{height}{tree height in meteres}
#'   \item{crownHeight}{crown height in meters}
#'   \item{BA}{basal area of individual trees in m2}
#'   \item{weight}{upscale weight to calculate hectare values}
#'   \item{p_BA}{basal area of individual trees in m2 from previous simulation step}
#'   \item{p_height}{tree height in meteres from previous simulation step}
#'   \item{p_crownHeight}{crown height in meters from previous simulation step}
#'   \item{p_weight}{upscale weight to calculate hectare values from previous simulation step}
#'   \item{volume}{tree volume in m3}
#'   \item{p_volume}{tree volume in m3 from previous simulation step}
#'
#' }
#' @export
"data_v3"
