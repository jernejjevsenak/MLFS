#' An example of data with individual tree and crown heights that can be used as
#' a fitting data frame for predicting tree and crown heights in MLFS
#'
#' This is simulated data that reassemble the national forest inventory  data.
#' We use it to show how to run examples for some specific functions
#'
#' @format A data frame with 2741 rows and 8 variables:
#' \describe{
#'   \item{plotID}{a unique identifier for plot}
#'   \item{treeID}{a unique identifier for tree}
#'   \item{year}{year in which plot was visited}
#'   \item{speciesGroup}{identifier for species group}
#'   \item{species}{species name}
#'   \item{height}{tree height in meters}
#'   \item{crownHeight}{crown height in meters}
#'   \item{BA}{basal area of individual trees in m2}
#' }
#'
#' @export
"data_tree_heights"
