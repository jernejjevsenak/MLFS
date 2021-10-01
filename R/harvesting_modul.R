#' simulate_harvesting
#'
#' Harvesting model for the MLFS
#' @keywords internal

simulate_harvesting <- function(df, harvesting_sum,
                                df_thinning_weights = NULL,
                                df_final_cut_weights = NULL,

                                df_thinning_weights_plot = NULL,
                                df_final_cut_weights_plot = NULL,

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
  volume_mid <- NULL
  weight <- NULL
  weight_mid <- NULL
  volume_ha <- NULL
  col_sum <- NULL
  code <- NULL
  plotID <- NULL
  treeID <- NULL
  protected <- NULL

  b <- df[df$code == 2,]
  a <- df[df$code != 2,]

  # merge thinning weights (species)
  if(!is.null(df_thinning_weights)){

    a <- merge(a, df_thinning_weights, by = 'species', all.x = TRUE)
    colnames(a)[ncol(a)] <- 'thinning_weight'

  } else {

    a$thinning_weight <- 1

  }

  # merge final_cut weights (species)
  if(!is.null(df_final_cut_weights)){

    a <- merge(a, df_final_cut_weights, by = 'species', all.x = TRUE)
    colnames(a)[ncol(a)] <- 'final_cut_weight'

  } else {

    a$final_cut_weight <- 1

  }





  # merge thinning weights (plot)
  if(!is.null(df_thinning_weights_plot)){

    a <- merge(a, df_thinning_weights_plot, by = 'plotID', all.x = TRUE)
    colnames(a)[ncol(a)] <- 'thinning_weight_plot'

  } else {

    a$thinning_weight_plot <- 1

  }

  # merge final_cut weights (plot)
  if(!is.null(df_final_cut_weights_plot)){

    a <- merge(a, df_final_cut_weights_plot, by = 'plotID', all.x = TRUE)
    colnames(a)[ncol(a)] <- 'final_cut_weight_plot'

  } else {

    a$final_cut_weight_plot <- 1

  }





if (harvest_sum_level == 1){  # regional (national level)

  if (plot_upscale_type == "factor"){

    a <- mutate(a, volume_ha = volume_mid * weight_mid * plot_upscale_factor)

  } else if (plot_upscale_type == "area") {

    n_plots <- length(unique(a$plotID))
    plot_represents <- forest_area_ha/n_plots

    a <- mutate(a, volume_ha = volume_mid * weight_mid * plot_represents)

  } else {

    stop("plot_upscale_type should be 'area' or 'factor'")

    }

  } else if (harvest_sum_level == 0){ # plot level

  a <- mutate(a, volume_ha = volume_mid * weight_mid)

  } else {

    stop("harvest_sum_level should be 0 or 1")

  }

  # 1 random harvesting
  if (harvesting_type == "random"){

    sampled_rows <- sample(1:NROW(a), size = nrow(a), prob = a$thinning_weight * a$thinning_weight_plot)
    a <- a[sampled_rows, ]
    a <- arrange(a, protected) # protected trees are at the bottom, so they won't be harvested

  } else if (harvesting_type == "final_cut"){

    sampled_rows <- sample(1:NROW(a), size = nrow(a), prob = (a$BA_mid ^ (final_cut_weight)) * a$final_cut_weight * a$final_cut_weight_plot)

    a <- a[sampled_rows, ]
    a <- arrange(a, protected) # protected trees are at the bottom, so they won't be harvested

    #################################################################


  } else if (harvesting_type == "thinning") {

    sampled_rows <- sample(1:NROW(a), size = nrow(a), prob = (a$BA_mid ^ (-thinning_small_weight)) * a$thinning_weight * a$thinning_weight_plot)
    a <- a[sampled_rows, ]
    a <- arrange(a, protected) # protected trees are at the bottom, so they won't be harvested

    ##############################################



  } else if (harvesting_type == "combined"){

    share_final_cut = 1 - share_thinning

    sum_volume <- sum(a$volume_mid, na.rm = T)

    aFC <- a
    aTH <- a

    #################
    # aFC Final_cut #
    #################

    sampled_rows <- sample(1:NROW(aFC), size = nrow(aFC), prob = (aFC$BA_mid ^ final_cut_weight) * aFC$final_cut_weight * aFC$final_cut_weight_plot)
    aFC <- aFC[sampled_rows, ]
    aFC <- arrange(aFC, protected)

    aFC <- mutate(aFC, col_sum = cumsum(replace_na(volume_ha, 0)),
                code = ifelse(col_sum < (harvesting_sum * share_final_cut), 1, code))

    aFC<- aFC[aFC$code == 1,]

    ################
    # aTH THINNING #
    ################

    #remove trees which were already cut
    aTH <- aTH[!(aTH$treeID %in% c(aFC$treeID)),]

    sampled_rows <- sample(1:NROW(aTH), size = nrow(aTH), prob = (aTH$BA_mid ^ -thinning_small_weight) * aTH$thinning_weight * aFC$thinning_weight_plot)
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

  a <- dplyr::select(a, colnames(df))

  df <- rbind(b, a)

  df <- arrange(df, plotID, treeID)

  return(df)
}




