#' calculate_standVars_halfPeriod
#'
#' Function for the calculation of stand variables: stand basal area and the
#' number of trees
#'
#' @keywords internal
#'

calculate_standVars_halfPeriod <- function(df){

  # Define global variables
  year <- NULL
  plotID <- NULL
  code <- NULL
  weight <- NULL
  BA <- NULL

  df$stand_BA_mid <- NULL
  df$stand_n_mid <- NULL

  data_stand <- group_by(df, year, plotID) %>%

    # dead/harvested/and ingrowth trees have reduced effect on stand variables
    mutate(weight = ifelse(code %in% c(1,2,3,15), weight /2, weight)) %>%

    summarise(stand_BA_mid = sum(BA_mid*weight_mid, na.rm = TRUE),
              stand_n_mid = sum(weight_mid, na.rm = T)
    )

  df <- merge(df, data_stand, by = c("year", "plotID"))

  return(df)

}
