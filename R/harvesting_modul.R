#' simulate_harvesting
#'
#' Harvesting model for the MLFS
#' @keywords internal

simulate_harvesting <- function(df, harvesting_sum,
                                harvesting_type = "random",
                                share_thinning = 0.8,
                                final_cut_weight = 10000000,
                                thinning_small_weight = 100000,

                                harvest_sum_level = 1,
                                plot_upscale_type,
                                plot_upscale_factor,
                                forest_area_ha

                                ){


  # Define global variables
  volume <- NULL
  weight <- NULL
  volume_ha <- NULL
  col_sum <- NULL
  code <- NULL
  plotID <- NULL
  treeID <- NULL
  protected <- NULL

  a <- df

if (harvest_sum_level == 1){  # regional (national level)

  if (plot_upscale_type == "factor"){

    a <- mutate(a, volume_ha = volume * weight * plot_upscale_factor)

  } else if (plot_upscale_type == "area") {

    n_plots <- length(unique(a$plotID))
    plot_represents <- forest_area_ha/n_plots

    a <- mutate(a, volume_ha = volume * weight * plot_represents)

  } else {

    stop("plot_upscale_type should be 'area' or 'factor'")

    }

  } else if (harvest_sum_level == 0){ # plot level

  a <- mutate(a, volume_ha = volume * weight)

  } else {

    stop("harvest_sum_level should be 0 or 1")

  }

  # 1 random harvesting
  if (harvesting_type == "random"){

    a <- a[sample(nrow(a)),]
    a <- arrange(a, protected) # protected trees are at the bottom, so they won't be harvested

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
    a <- arrange(a, protected) # protected trees are at the bottom, so they won't be harvested

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
    a <- arrange(a, protected) # protected trees are at the bottom, so they won't be harvested

  } else if (harvesting_type == "combined"){

    share_final_cut = 1 - share_thinning

    sum_volume <- sum(a$volume, na.rm = T)

    aFC <- a
    aTH <- a

    #################
    # aFC Final_cut #
    #################

    aFC$share_volume <- aFC$volume / sum_volume

    # more intense
    aFC$share_volume <- ifelse(aFC$share_volume > mean(aFC$share_volume, na.rm = TRUE),
                               aFC$share_volume * final_cut_weight,
                               aFC$share_volume * (1/final_cut_weight))

    sampled_rows <- sample(1:NROW(aFC), size = nrow(aFC), prob = aFC$share_volume)
    aFC <- aFC[sampled_rows, ]
    aFC <- arrange(aFC, protected)

    aFC <- mutate(aFC, col_sum = cumsum(replace_na(volume_ha, 0)),
                code = ifelse(col_sum < (harvesting_sum * share_final_cut), 1, code))

    aFC<- aFC[aFC$code == 1,]

    ################
    # aTH THINNING #
    ################

    # aFC_cut <- aFC[aFC$code == 1, "treeID"]

    #remove trees which were already cut
    aTH <- aTH[!(aTH$treeID %in% c(aFC$treeID)),]
    aTH$share_volume <- (1 - (aTH$volume / sum_volume))

    # more intense
    aTH$share_volume <- ifelse(aTH$share_volume > mean(aTH$share_volume, na.rm = TRUE),
                             aTH$share_volume * thinning_small_weight,
                             aTH$share_volume * (1/thinning_small_weight))

    sampled_rows <- sample(1:NROW(aTH), size = nrow(aTH), prob = aTH$share_volume)
    aTH <- mutate(aTH, col_sum = cumsum(replace_na(volume_ha, 0)),
                  code = ifelse(col_sum < (harvesting_sum * share_thinning), 1, code))

    aTH <- aTH[sampled_rows, ]
    aTH <- arrange(aTH, protected) # protected trees are at the bottom, so they won't be harvested

    #############

    a <- rbind(aTH, aFC)


  } else {

    stop(paste0("harvesting_type should be one of 'random', 'final_cut' or 'thinning', but instead it is ", harvesting_type ))

  }

  if (harvesting_type != "combined"){

    a <- mutate(a, col_sum = cumsum(replace_na(volume_ha, 0)),
                code = ifelse(col_sum < harvesting_sum, 1, code))

  }

  df <- select(a, colnames(df))
  df <- arrange(df, plotID, treeID)

  return(df)
}




