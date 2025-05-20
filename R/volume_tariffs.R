#' volume_tariffs
#'
#' One-parameter volume functions (tariffs) for the MLFS.
#'
#' @param df data frame with tree heights and basal areas for individual trees
#' @param data_tariffs data frame with plot- and species-specific parameters for
#' the calculations of tree volume
#'
#' @return a data frame with calculated volume for all trees
#'
#' @examples
#' data(data_v3)
#' data(data_tariffs)
#' data_v3 <- volume_tariffs(df = data_v3, data_tariffs = data_tariffs)

volume_tariffs <- function(df, data_tariffs){

# Define global variables
BA <- NULL
p_BA <- NULL
tarifa_class <- NULL

  initial_colnames <- colnames(df)

    # unique(df$species[!  (df$species %in% data_tariffs$species)])

  df <- merge(df, data_tariffs, by = c("plotID", "species", "year"), all.x = TRUE)

  test_df <- dplyr::filter(df, is.na(v45)) %>% dplyr::select(plotID, species, year)

  if (nrow(test_df) > 0) {

    # 1) catch any actual NA’s in species or plotID
    na_rows <- which(is.na(test_df$species) | is.na(test_df$plotID))
    if (length(na_rows) > 0) {
      msgs1 <- sapply(na_rows, function(i) {
        sp  <- ifelse(is.na(test_df$species[i]), "<MISSING_SPECIES>", test_df$species[i])
        pid <- ifelse(is.na(test_df$plotID[i]),   "<MISSING_PLOTID>",   test_df$plotID[i])
        paste0("row ", i, ": species='", sp, "' @ plotID='", pid, "'")
      })
      stop(
        paste0(
          "Missing data in data_tariffs (NA’s): ",
          paste(msgs1, collapse = "; ")
        )
      )
    }

    # 2) build the full grid of all observed species × plotID
    full_grid <- expand.grid(
      plotID  = unique(test_df$plotID),
      species = unique(test_df$species),
      stringsAsFactors = FALSE
    )

    # pull out just the species/plot pairs you actually have
    observed <- unique(test_df[, c("plotID", "species")])

    # find the combos in full_grid that never occur in observed
    missing_grid <- merge(full_grid, observed,
                          by = c("plotID", "species"),
                          all.x = TRUE
    )
    missing_grid <- missing_grid[ ! (paste0(missing_grid$plotID, "|", missing_grid$species)
                                     %in% paste0(observed$plotID, "|", observed$species)),
    ]

    if (nrow(missing_grid) > 0) {
      msgs2 <- with(missing_grid,
                    paste0("species='", species, "' @ plotID='", plotID, "'"))
      stop(
        paste0(
          "Missing data in data_tariffs (absent combos): ",
          paste(msgs2, collapse = "; ")
        )
      )
    }
  }

  # df_ <- dplyr::filter(df_, is.na(v45))
  # write.csv(df_, "missing_tariffs.csv", row.names = FALSE)

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
