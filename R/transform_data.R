#' transform_data
#'
#' Initial data transformation
#'
#' @keywords internal

transform_data <- function(df, include_climate, df_climate, select_months_climate){

  # Define global variables
  year <- NULL
  plotID<- NULL
  treeID<- NULL
  BA<- NULL
  p_BA<- NULL
  month<- NULL
  p_sum<- NULL
  t_avg<- NULL
  height <- NULL
  crownHeight <- NULL

  # Next step: create DF for yield and predict heights and calculate volume
  unique_years <- sort(unique(df$year), decreasing = T)

  listed <- list()
  b = 1

  for (i in 1:length(unique_years)){

    if (is.na(unique_years[i + 1]))
      next()

    df_temp <- dplyr::filter(df, year == unique_years[i])
    df_temp_year_before <- dplyr::filter(df, year == unique_years[i+1])
    df_temp_year_before <- select(df_temp_year_before, plotID, treeID, BA, height, crownHeight)
    colnames(df_temp_year_before)[3:5] <- c("p_BA", "p_height", "p_crownHeight")
    df_temp <- merge(df_temp, df_temp_year_before, by = c("plotID", "treeID"), all.x = TRUE)
    df_temp <- mutate(df_temp, BAI = BA - p_BA)

    if (include_climate == TRUE){

      max_year <- max(df_temp$year)
      min_year <- unique_years[i + 1]

      climate_fit <- filter(df_climate, year %in% seq(min_year, max_year)) %>% filter(month %in% select_months_climate) %>%
        group_by(plotID) %>% summarise(p_sum = sum(p_sum), t_avg = mean(t_avg))

      df_temp <- merge(df_temp, climate_fit, by = "plotID")

    } else {

      df_temp$p_sum <- NA
      df_temp$t_avg <- NA

    }

    listed[[b]] <- df_temp
    b = b + 1

  }

  df <- do.call(rbind, listed)

  return(df)

}
