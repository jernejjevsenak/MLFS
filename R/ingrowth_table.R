#' An example data of ingrowth_table
#'
#' Ingrowth table is used within the ingrowth sub model to correctly simulate
#' different ingrowth levels and associated upscale weights
#'
#' @format A data frame with 2 rows and 4 variables:
#' \describe{
#'   \item{code}{ingrowth codes}
#'   \item{DBH_threshold}{a DBH threshold for particular ingrowth category}
#'   \item{DBH_max}{maximum DBH for a particular ingrowth category}
#'   \item{weight}{the upscale weight for particular measurement category}
#' }
#' @export
"ingrowth_table"
