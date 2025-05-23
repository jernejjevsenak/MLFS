#' An example of joined national forest inventory and site data that is used
#' within the MLFS
#'
#' This is simulated data that reassemble the national forest inventory data.
#' We use it to show how to run examples for simulating Basal Area Increments
#' (BAI) and the ingrowth of new trees. To make examples running more quickly,
#' we keep only one tree species: PINI
#'
#' @format A data frame with 186 rows and 27 variables:
#' \describe{
#'   \item{species}{species name}
#'   \item{year}{year in which plot was visited}
#'   \item{plotID}{a unique identifier for plot}
#'   \item{treeID}{a unique identifier for tree}
#'   \item{speciesGroup}{identifier for species group}
#'   \item{code}{status of a tree: 0 (normal), 1(harvested), 2(dead), 3 (ingrowth)}
#'   \item{height}{tree height in meters}
#'   \item{crownHeight}{crown height in meters}
#'   \item{protected}{logical, 1 if protected, otherwise 0}
#'   \item{slope}{slope on a plot}
#'   \item{elevation}{plot elevation}
#'   \item{northness}{plot northness, 1 is north, 0 is south}
#'   \item{siteIndex}{a proxy for site index, higher value represents more productive sites}
#'   \item{BA}{basal area of individual trees in m2}
#'   \item{weight}{upscale weight to calculate hectare values}
#'   \item{stand_BA}{Total stand basal area}
#'   \item{stand_n}{The number of trees in a stand}
#'   \item{BAL}{Basal Area in Large trees}
#'   \item{p_BA}{basal area of individual trees in m2 from previous simulation step}
#'   \item{p_height}{tree height in meters from previous simulation step}
#'   \item{p_crownHeight}{crown height in meters from previous simulation step}
#'   \item{p_weight}{upscale weight to calculate hectare values from previous simulation step}
#'   \item{BAI}{basal area increment}
#'   \item{p_sum}{monthly precipitation sum}
#'   \item{t_avg}{monthly mean temperature}
#'   \item{volume}{tree volume in m3}
#'   \item{p_volume}{tree volume in m3 from previous simulation step}
#' }
#' @export
"data_v6"
