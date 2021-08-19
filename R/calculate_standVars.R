#' calculate_standVars
#'
#' Function for the calculation of stand variables: stand basal area and the
#' number of trees
#'
#' @keywords internal
#'

calculate_standVars <- function(df){

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
