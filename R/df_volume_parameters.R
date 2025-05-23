#' An example table with parameters and equations for n-parametric volume functions
#'
#' Volume functions can be specified for each species and plot separately, also limited to specific DBH interval. The factor variables (vol_factor,
#' h_factor and DBH_factor) are used to control the input and output units.
#'
#' @format A data frame with 6 rows and 14 variables:
#' \describe{
#'   \item{species}{species name as used in data_NFI. The category REST is used for all species without specific equation}
#'   \item{equation}{equation for selected volume function}
#'   \item{vol_factor}{will be multiplied with the volume}
#'   \item{h_factor}{will be multiplied with tree height}
#'   \item{d_factor}{will be divided with tree DBH}
#'   \item{DBH_min}{lower interval threshold for considered trees}
#'   \item{DBH_max}{upper interval threshold for considered trees}
#'   \item{a}{parameter a for volume equation}
#'   \item{b}{parameter b for volume equation}
#'   \item{c}{parameter c for volume equation}
#'   \item{d}{parameter d for volume equation}
#'   \item{e}{parameter e for volume equation}
#'   \item{f}{parameter f for volume equation}
#'   \item{g}{parameter g for volume equation}
#' }
#' @export
"df_volume_parameters"
