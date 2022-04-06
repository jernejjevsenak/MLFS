#' An example of data_ingrowth suitable for the MLFS
#'
#' An example of plot-level data with plotID, stand variables and site
#' descriptors, and the two target variables describing the number of ingrowth
#' trees for inner (ingrowth_3) and outer (ingrowth_15) circles
#'
#' @format A data frame with 365 rows and 11 variables:
#' \describe{
#'   \item{plotID}{a unique identifier for plot}
#'   \item{year}{year in which plot was visited}
#'   \item{stand_BA}{Total stand basal area}
#'   \item{stand_n}{The number of trees in a stand}
#'   \item{BAL}{Basal Area in Large trees}
#'   \item{slope}{slope on a plot}
#'   \item{elevation}{plot elevation}
#'   \item{siteIndex}{a proxy for site index, higher value represents more productible sites}
#'   \item{northness}{plot northness, 1 is north, 0 is south}
#'   \item{ingrowth_3}{the number of new trees in inner circle}
#'   \item{ingrowth_15}{the number of new trees in outer circle}
#' }
#' @export
"data_ingrowth"
