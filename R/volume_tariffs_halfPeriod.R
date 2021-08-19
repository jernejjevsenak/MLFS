#' vol_tariffs_halfPeriod
#'
#' One-parameter volume functions (tariffs) for the MLFS (half period)
#' @keywords internal

vol_tariffs_halfPeriod <- function(df, data_tariffs){

# Define global variables
BA_mid <- NULL
tarifa_class <- NULL

  initial_colnames <- colnames(df)

  df <- merge(df, data_tariffs, by = c("plotID", "species"))

  df <- mutate(df,
               D = sqrt(4*BA_mid/pi) * 100,
               tarifa_class = as.numeric(tarifa_class)
  )

  df$volume_mid <- ifelse(df$tarifa_class <= 20 & df$D >= 25.0, df$v45 / 1400 * (df$D - 5) * (df$D - 10),
                      ifelse(df$tarifa_class <= 20 & df$D < 25.0, df$v45 / 1400.0 * (-226.33 + 38.575*df$D - 1.9237 * (df$D)^2 + 0.04876 * (df$D)^3),
                             ifelse(df$tarifa_class > 20 & df$tarifa_class <= 40, df$v45 / 1600.0 * (df$D - 2.5) * (df$D - 7.5),
                                    df$v45 / 1800 * df$D * (df$D - 5))))

  df <- select(df, all_of(initial_colnames))

  return(df)
}
