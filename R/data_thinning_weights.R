#' An example of data_thinning_weights
#'
#' Each species should have one weight that is multiplied with the probability of being harvested when thinning is applied
#'
#' @format A data frame with 36 rows and 6 variables:
#' \describe{
#'   \item{species}{species name as used in data_NFI}
#'   \item{step_1}{thinning weight applied in step 1}
#'   \item{step_2}{thinning weight applied in step 2}
#'   \item{step_3}{thinning weight applied in step 3}
#'   \item{step_4}{thinning weight applied in step 4}
#'   \item{step_5}{thinning weight applied in step 5 and all subsequent steps}
#' }
#' @export
"data_thinning_weights"
