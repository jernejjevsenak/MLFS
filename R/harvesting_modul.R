#' A sub model to simulate harvesting within the MLFS
#'
#' Harvesting is based on probability sampling, which depends on the selected
#' parameters and the seize of a tree. Bigger trees have higher probability of
#' being harvested when final cut is applied, while smaller trees have higher
#' probability of being sampled in the case of thinning.
#'
#' @param df a data frame with individual tree data, which include basal areas
#' in the middle of a simulation step, species name and code
#' @param harvesting_sum a value, or a vector of values defining the harvesting
#' sums through the simulation stage. If a single value, then it is used in all
#' simulation steps. If a vector of values, the first value is used in the first
#' step, the second in the second step, etc.
#' @param df_thinning_weights_species data frame with thinning weights for each
#' species. The first column represents species code, each next column consists
#' of species-specific thinning weights
#' @param df_final_cut_weights_species data frame with final cut weights for each
#' species. The first column represents species code, each next column consists
#' of species-specific final cut weights
#' @param df_thinning_weights_plot data frame with harvesting weights related
#' to plot IDs, used for thinning
#' @param df_final_cut_weights_plot data frame with harvesting weights related
#' to plot IDs, used for final cut
#' @param harvesting_type character, it could be 'random', 'final_cut',
#' 'thinning' or 'combined'. The latter combines 'final_cut' and 'thinning'
#' options, where the share of each is specified with the argument
#' 'share_thinning'
#' @param share_thinning numeric, a number between 0 and 1 that specifies the
#' share of thinning in comparison to final_cut. Only used if harvesting_type
#' is 'combined'
#' @param final_cut_weight numeric value affecting the probability distribution
#' of harvested trees. Greater value increases the share of harvested trees
#' having larger DBH. Default is 10.
#' @param thinning_small_weight numeric value affecting the probability
#' distribution of harvested trees. Greater value increases the share of
#' harvested trees having smaller DBH. Default is 1.
#' @param harvest_sum_level integer with value 0 or 1 defining the level of
#' specified harvesting sum: 0 for plot level and 1 for regional level
#' @param plot_upscale_type character defining the upscale method of plot level
#' values. It can be 'area' or 'upscale factor'. If 'area', provide the forest
#' area represented by all plots in hectares (forest_area_ha argument). If
#' 'factor', provide the fixed factor to upscale the area of all plots. Please
#' note: forest_area_ha/plot_upscale_factor = number of unique plots. This
#' argument is important when harvesting sum is defined on regional level.
#' @param plot_upscale_factor numeric value to be used to upscale area of each
#' plot
#' @param forest_area_ha the total area of all forest which are subject of the
#' simulation
#'
#' @return a data frame with updated status (code) of all individual trees based
#' on the simulation of harvesting
#'
#' @examples
#'
#' library(MLFS)
#' data(data_v5)
#'
#' data_v5 <- simulate_harvesting(df = data_v5,
#'             harvesting_sum = 5500000,
#'             harvesting_type = "combined",
#'             share_thinning = 0.50,
#'             harvest_sum_level = 1,
#'             plot_upscale_type = "factor",
#'             plot_upscale_factor = 1600,
#'             final_cut_weight = 5,
#'             thinning_small_weight = 1)
#'

simulate_harvesting <- function(df, harvesting_sum,

                                df_thinning_weights_species = NULL,
                                df_final_cut_weights_species = NULL,

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
  if(!is.null(df_thinning_weights_species)){

    a <- merge(a, df_thinning_weights_species, by = 'species', all.x = TRUE)
    colnames(a)[ncol(a)] <- 'thinning_weight_species'

    a$thinning_weight_species <- ifelse(is.na(a$thinning_weight_species), 1, a$thinning_weight_species)

    # Probabilities can not be 0
    a$thinning_weight_species <- ifelse(a$thinning_weight_species < 0.00001, 0.001, a$thinning_weight_species)

  } else {

    a$thinning_weight_species <- 1

  }

  # merge final_cut weights (species)
  if(!is.null(df_final_cut_weights_species)){

    a <- merge(a, df_final_cut_weights_species, by = 'species', all.x = TRUE)
    colnames(a)[ncol(a)] <- 'final_cut_weight_species'

    a$final_cut_weight_species <- ifelse(is.na(a$final_cut_weight_species), 1, a$final_cut_weight_species)
    # Probabilities can not be 0
    a$final_cut_weight_species <- ifelse(a$final_cut_weight_species < 0.00001, 0.001, a$final_cut_weight_species)


  } else {

    a$final_cut_weight_species <- 1

  }


  # merge thinning weights (plot)
  if(!is.null(df_thinning_weights_plot)){

    a <- merge(a, df_thinning_weights_plot, by = 'plotID', all.x = TRUE)
    colnames(a)[ncol(a)] <- 'thinning_weight_plot'

    a$thinning_weight_plot <- ifelse(is.na(a$thinning_weight_plot), 1, a$thinning_weight_plot)

    # Probabilities can not be 0
    a$thinning_weight_plot <- ifelse(a$thinning_weight_plot < 0.00001, 0.001, a$thinning_weight_plot)

  } else {

    a$thinning_weight_plot <- 1

  }

  # merge final_cut weights (plot)
  if(!is.null(df_final_cut_weights_plot)){

    a <- merge(a, df_final_cut_weights_plot, by = 'plotID', all.x = TRUE)
    colnames(a)[ncol(a)] <- 'final_cut_weight_plot'

    a$final_cut_weight_plot <- ifelse(is.na(a$final_cut_weight_plot), 1, a$final_cut_weight_plot)
    # Probabilities can not be 0
    a$final_cut_weight_plot <- ifelse(a$final_cut_weight_plot < 0.00001, 0.001, a$final_cut_weight_plot)

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

    sampled_rows <- sample(1:NROW(a), size = nrow(a), prob = 1 * a$thinning_weight_species * a$thinning_weight_plot)
    a <- a[sampled_rows, ]
    a <- arrange(a, protected) # protected trees are at the bottom, so they won't be harvested

  } else if (harvesting_type == "final_cut"){

    sampled_rows <- sample(1:NROW(a), size = nrow(a), prob = (a$BA_mid ^ (final_cut_weight)) * a$final_cut_weight_species * a$final_cut_weight_plot)

    a <- a[sampled_rows, ]
    a <- arrange(a, protected) # protected trees are at the bottom, so they won't be harvested

    #################################################################

  } else if (harvesting_type == "thinning") {

    sampled_rows <- sample(1:NROW(a), size = nrow(a), prob = (a$BA_mid ^ (-thinning_small_weight)) * a$thinning_weight_species * a$thinning_weight_plot)
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

    sampled_rows <- sample(1:NROW(aFC), size = nrow(aFC), prob = (aFC$BA_mid ^ final_cut_weight) * aFC$final_cut_weight_species * aFC$final_cut_weight_plot)
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

    sampled_rows <- sample(1:NROW(aTH), size = nrow(aTH), prob = (aTH$BA_mid ^ -thinning_small_weight) * aTH$thinning_weight_species * aTH$thinning_weight_plot)

    aTH <- aTH[sampled_rows, ]
    aTH <- arrange(aTH, protected) # protected trees are at the bottom, so they won't be harvested

    aTH <- mutate(aTH, col_sum = cumsum(replace_na(volume_ha, 0)),
                  code = ifelse(col_sum < (harvesting_sum * share_thinning), 1, code))

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




