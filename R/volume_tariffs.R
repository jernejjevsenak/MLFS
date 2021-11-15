#' vol_tariffs
#'
#' One-parameter volume functions (tariffs) for the MLFS
#' @keywords internal

vol_tariffs <- function(df, data_tariffs){

# Define global variables
BA <- NULL
p_BA <- NULL
tarifa_class <- NULL

  initial_colnames <- colnames(df)

  df <- merge(df, data_tariffs, by = c("plotID", "species"))

  df <- mutate(df,
               D = sqrt(4*BA/pi) * 100,
               p_D = sqrt(4*p_BA/pi) * 100,
               tarifa_class = as.numeric(tarifa_class)
  )

  df$volume <- ifelse(df$tarifa_class <= 20 & df$D >= 25.0, df$v45 / 1400 * (df$D - 5) * (df$D - 10),
                      ifelse(df$tarifa_class <= 20 & df$D < 25.0, df$v45 / 1400.0 * (-226.33 + 38.575*df$D - 1.9237 * (df$D)^2 + 0.04876 * (df$D)^3),
                             ifelse(df$tarifa_class > 20 & df$tarifa_class <= 40, df$v45 / 1600.0 * (df$D - 2.5) * (df$D - 7.5),
                                    df$v45 / 1800 * df$D * (df$D - 5))))

  df$p_volume <- ifelse(df$tarifa_class <= 20 & df$p_D >= 25.0, df$v45 / 1400 * (df$p_D - 5) * (df$p_D - 10),
                        ifelse(df$tarifa_class <= 20 & df$p_D < 25.0, df$v45 / 1400.0 * (-226.33 + 38.575*df$p_D - 1.9237 * (df$p_D)^2 + 0.04876 * (df$p_D)^3),
                               ifelse(df$tarifa_class > 20 & df$tarifa_class <= 40, df$v45 / 1600.0 * (df$p_D - 2.5) * (df$p_D - 7.5),
                                      df$v45 / 1800 * df$p_D * (df$p_D - 5))))

  df <- dplyr::select(df, all_of(initial_colnames))

  return(df)
}
