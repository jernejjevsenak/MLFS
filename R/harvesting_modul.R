#' simulate_harvesting
#'
#' Harvesting model for the MLFS
#' @keywords internal

simulate_harvesting <- function(df, harvesting_sum, forest_area_ha,
                                harvesting_type = "random",
                                final_cut_weight = 1, thinning_small_weight = 10){


# Define global variables
volume <- NULL
weight <- NULL
volume_ha <- NULL
col_sum <- NULL
code <- NULL
plotID <- NULL
treeID <- NULL

  a <- df

  n_plots <- length(unique(a$plotID))
  plot_represents <- forest_area_ha/n_plots

  a <- mutate(a, volume_ha = volume * weight * plot_represents)

  # 1 random harvesting
  if (harvesting_type == "random"){

    a <- a[sample(nrow(a)),]

  } else if (harvesting_type == "final_cut"){

    # 1 final cut
    sum_volume <- sum(a$volume, na.rm = T)
    a$share_volume <- a$volume / sum_volume

    # more intense
    a$share_volume <- ifelse(a$share_volume > mean(a$share_volume, na.rm = TRUE),
                             a$share_volume * final_cut_weight,
                             a$share_volume * (1/final_cut_weight))

    sampled_rows <- sample(1:NROW(a), size = nrow(a), prob = a$share_volume)

    a <- a[sampled_rows, ]

  } else if (harvesting_type == "thinning") {

    # 3 thinning which prefers smaller trees
    sum_volume <- sum(a$volume, na.rm = T)
    a$share_volume <- (1 - (a$volume / sum_volume))

    # more intense
    a$share_volume <- ifelse(a$share_volume > mean(a$share_volume, na.rm = TRUE),
                             a$share_volume * thinning_small_weight,
                             a$share_volume * (1/thinning_small_weight))

    sampled_rows <- sample(1:NROW(a), size = nrow(a), prob = a$share_volume)

    a <- a[sampled_rows, ]

  } else {

    stop(paste0("harvesting_type should be on of 'random', 'final_cut' or 'thinning', but it is ", harvesting_type ))

  }

  a <- mutate(a, col_sum = cumsum(replace_na(volume_ha, 0)),
              code = ifelse(col_sum < harvesting_sum, 1, code))

  df <- select(a, colnames(df))
  df <- arrange(df, plotID, treeID)

  return(df)
}




