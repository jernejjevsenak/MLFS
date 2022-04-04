#' add_stand_variables
#'
#' This function adds two variables to existing data frame of individual tree
#' measurements: 1) stand basal area and 2) the number of trees per hectare
#'
#' @param df a data frame with individual tree measurements that include basal
#' area and the upscale factors. All trees should also be described with plotID
#' and year variables
#'
#' @examples
#' data(data_v1)
#' data_v1 <- add_stand_variables(df = data_v1)
#'

add_stand_variables <- function(df){

  # Define global variables
  year <- NULL
  plotID <- NULL
  code <- NULL
  weight <- NULL
  BA <- NULL

  df$stand_BA <- NULL
  df$stand_n <- NULL

  data_stand <- group_by(df, year, plotID) %>%

    # harvested trees have reduced effect on stand variables
    mutate(weight = ifelse(code %in% c(1), weight /2, weight)) %>%

    summarise(stand_BA = sum(BA*weight, na.rm = TRUE),
              stand_n = sum(weight, na.rm = T)
    )

  df <- merge(df, data_stand, by = c("year", "plotID"))

  return(df)

}
