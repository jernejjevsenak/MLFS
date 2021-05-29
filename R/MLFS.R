#' MLFS
#'
#' Machine Learning Forest Simulator
#' @param data_NFI data frame with individual tree variables
#' @param data_site data frame with site descriptors. This data is related to
#' data_NFI based on the 'plotID' column
#' @param data_tariffs optional, but mandatory if volume is calculated using the
#' one-parametric tariff functions. Data frame with plotID, species and V45. See
#' details.
#' @param data_climate data frame with climate data, covering the initial
#' calibration period and all the years which will be included in the simulation
#' @param sim_mortality logical, should mortality be simulated?
#' @param data_volF_param optional, data frame with species-specific volume
#' function parameters
#' @param form_factors optional, data frame with species-specific form factors
#' @param volume_calculation character string defining the method for volume
#' calculation: 'tariffs', 'volume_functions', 'form_factors' or
#' 'slo_2p_volume_functions'
#' @param sim_harvesting logical, should harvesting be simulated?
#' @param harvest_sum_level integer with value 0 or 1 defining the level of
#' specified harvesting sum: 0 for plot level and 1 for regional level.
#' @param plot_upscale_type character defining the upscale method of plot level
#' area. It can be 'area' or 'upscale factor'. If 'area', provide the forest
#' area represented by all plots in hectars (forest_area_ha argument). If
#' 'factor', provide the fixed factor to upscale the area of all plots. Please
#' note: forest_area_ha/plot_upscale_factor = number of unique plots. This
#' argument is important when harvesting sum is defined on regional level.
#' @param plot_upscale_factor numeric value to be used to upscale area of each
#' plot
#' @param harvesting_sum a value, or a vector of values defining the harvesting
#' sums through the simulation stage. If a single value, then it is used in all
#' simulation steps. If a vector of values, the first value is used in the first
#' step, the second in the second step, etc.
#' @sim_mortality logical, should mortality be simulated?
#' @param mortality_share a value, or a vector of values defining the proportion
#' of the volume which is to be the subject of mortality. If a single value,
#' then it is used in all simulation steps. If a vector of values, the first
#' value is used in the first step, the second in the second step, and so on.
#' @param forest_area_ha the total area of all forest which are subject of the
#' simulation.
#' @param harvesting_type character, it could be 'random', 'final cut' or
#' 'thinning'
#' @param final_cut_weight numeric value affecting the probability distribution
#' of harvested trees. Greater value increases the share of harvested trees
#' having larger DBH. Default is 1.
#' @param thinning_small_weight numeric value affecting the probability
#' distribution of harvested trees. Greater value increases the share of
#' harvested trees having smaller DBH. Default is 1.
#' @param k the number of folds to be used in the k fold cross-validation
#' @param blocked_cv logical, should the blocked cross-validation be used
#' @param species_n_threshold a positive integer defining the minimum number of
#' observations required to treat a species as an independent group
#' @param height_model character string defining the model to be used for height
#' prediction. If brnn, then ANN method with Bayesian Regularization is applied.
#' In addition, all 2- and 3- parametric H-D models from lmfor R package are
#' available.
#' @param crownHeight_model character string defining the model to be used for
#' crown heights. Available are ANN with Bayesian regularization (brnn) or
#' linear regression (lm)
#' @param BRNN_neurons positive integer defining the number of neurons to bo
#' used in the brnn method.
#' @param mortality_model model to be used for mortality prediction: 'glm' for
#' generalized linear models; 'rf' for random forest algorithm; 'naiveBayes' for
#' Naive Bayes algorithm
#' @param nb_laplace value used for Laplace smoothing (additive smoothing) in
#' naive Bayes algorithm. Defaults to 0 (no Laplace smoothing).
#' @param rf_mtry number of variables randomly sampled as candidates at each
#' split of a random forest model. If NULL, default settings are applied.
#' @param ingrowth_model model to be used for ingrowth predictions. 'glm' for
#' generalized linear models and 'rf' for random forest
#' @param sim_steps The number of simulation steps
#' @param merchantable_whole_tree character, 'merchantable' or 'whole_tree'. It
#' indicates which type of volume functions will be used. This parameter is used
#' only for volume calculation using the 'slo_2p_volume_functions'.
#' @param height_pred_level integer with value 0 or 1 defining the level of
#' prediction for height-diameter (H-D) models. The value 1 defines a plot-level
#' prediction, while the value 0 defines regional-level predictions. Default is
#' 0. If using 1, make sure to have representative plot-level data for each
#' species.
#' @param include_climate logical, should climate variables be included as
#' predictors
#' @param select_months_climate vector of subset months to be considered.
#' Default is c(1:12), which uses all months.
#' @param set_eval_mortality logical, should the mortality model be evaluated
#' and returned as the output
#' @param set_eval_crownHeight logical, should the crownHeight model be
#' evaluated and returned as the output
#' @param set_eval_height logical, should the height model be evaluated and
#' returned as the output
#' @param set_eval_ingrowth logical, should the the ingrowth model be evaluated
#' and returned as the output
#' @param set_eval_BAI logical, should the the BAI model be evaluated and
#' returned as the output

MLFS <- function(data_NFI, data_site,
                 data_tariffs = NULL,
                 data_climate = NULL,
                 data_volF_param = NULL,
                 form_factors = NULL, sim_steps,
                 volume_calculation = "volume_functions",
                 merchantable_whole_tree = "merchantable",
                 sim_harvesting = TRUE, sim_mortality = TRUE,
                 harvesting_sum,

                 forest_area_ha,
                 harvest_sum_level = 1,
                 plot_upscale_type,
                 plot_upscale_factor,

                 mortality_share = NA,
                 mortality_model = "rf",
                 ingrowth_model = "glm",
                 rf_mtry = NULL,
                 nb_laplace = 0,
                 harvesting_type = "final_cut", final_cut_weight = 1,
                 species_n_threshold = 100,
                 thinning_small_weight = 10,
                 height_model = "naslund",
                 crownHeight_model = "lm",
                 BRNN_neurons = 3,
                 height_pred_level = 0, # prediction level for lmfor (0 regional, 1 plot level) - works only if you have all species on all plots
                 include_climate = FALSE, select_months_climate = c(1:12),
                 set_eval_mortality = TRUE,
                 set_eval_crownHeight = TRUE,
                 set_eval_height = TRUE,
                 set_eval_ingrowth = TRUE,
                 set_eval_BAI = TRUE,
                 k = 10, blocked_cv = TRUE){

  # Define global variables
  DBH <- NULL
  height <- NULL
  crownHeight <- NULL
  BAI <- NULL
  BA <- NULL
  code <- NULL
  plotID <- NULL
  year <- NULL
  stand_BA <- NULL
  stand_n <- NULL
  BAL <- NULL
  ingrowth_small <- NULL
  ingrowth_big <- NULL

  pb = txtProgressBar(min = 0, max = sim_steps, initial = 0, style = 3)

  setTxtProgressBar(pb, 1)

  options(dplyr.summarise.inform= FALSE)

  # sim_steps should be positive
  if (sim_steps < 1){
    stop("sim_steps should be at least 1")
  }

  # What is the simulation step
  # I am also considering more than 2 calibration periods - we take the last two
  sim_step_years <- sort(unique(as.numeric(data_NFI$year)))
  sim_step_years <- sim_step_years[length(sim_step_years)] - sim_step_years[length(sim_step_years)-1]

  # check the first column
  if (colnames(data_site)[1] != "plotID"){

    stop(paste0("The first column of data_site should be 'plotID', but instead it is ", colnames(data_site)[1]))

  }

  # save site variable names and use them in formulas
  site_vars <- colnames(data_site)[-1] # plotID should be removed




  # merge NFI and site descriptors
  data <- merge(data_NFI, data_site, by = "plotID")

  # 1 Calculate stand basal area and number of trees
  data <- dplyr::mutate(data,
                        BA = ((DBH/2)^2 * pi)/10000,
                        DBH = NULL
                        )

  # kode iz NFI so: 0 (normal), 1 (harvested),  2 (dead), 3 (ingrowth), 15 (ingrowth)?
  data <- calculate_standVars(df = data)

  # specify form factors
  if (is.null(form_factors)){

    data$form <- 0.42

  } else {

    data <- merge(data, form_factors, by = "species")

  } # You can also add here third option, merge by plot & species

  # 2 Calculate BAL
  data$BAL <- NA
  data <- calculate_BAL(df = data)

  # 3 Transform data, where yield is properly expressed
  data <- transform_data(df = data, include_climate = include_climate,
                         df_climate = data_climate, select_months_climate = select_months_climate)

  ######################################################
  # create fitting data frames for
  # 1) height
  data_height <- dplyr::filter(data, !is.na(height))

  # 2) crown height
  data_crownHeight <- dplyr::filter(data, !is.na(crownHeight))

  # 3) BAI
  data_BAI <- dplyr::filter(data, !is.na(BAI))

  # 4) mortality
  data$p_height <- NA
  data$p_crownHeight <- NA
  data_mortality <- data

  h_predictions <- height_prediction(df_fit = data_mortality,
                                      df_predict = data_mortality,
                                      species_n_threshold = species_n_threshold,
                                      height_pred_level = height_pred_level,
                                      height_model = height_model,
                                      BRNN_neurons = BRNN_neurons,
                                      eval_model_height = set_eval_height,
                                      blocked_cv = blocked_cv, k = k)

  data_mortality <- h_predictions$data_height_predictions

  Crown_h_predictions <- crownHeight_prediction(df_fit = data_mortality,
                                           df_predict = data_mortality,
                                           site_vars = site_vars,
                                           crownHeight_model = crownHeight_model,
                                           BRNN_neurons = BRNN_neurons,
                                           species_n_threshold = species_n_threshold,
                                           k = k, blocked_cv = blocked_cv,
                                           eval_model_crownHeight = set_eval_crownHeight)

  data_mortality <- Crown_h_predictions$predicted_crownHeight

  data_mortality <- dplyr::filter(data_mortality, !is.na(BA))

  # 5) Ingrowth - mora biti brez samo kode 0 3 15
  data_ingrowth <- data_mortality
  data_ingrowth <- dplyr::filter(data_ingrowth, code %in% c(0, 3, 15)) %>%
    mutate(ingrowth_small = ifelse(code == 3, 1, 0),
           ingrowth_big = ifelse(code == 15, 1, 0)
    )

  data_ingrowth_stand <- select(data_ingrowth, plotID, year, stand_BA, stand_n, BAL, all_of(site_vars)) %>%
    group_by(plotID, year) %>% summarise_all(.funs = mean, na.rm = TRUE)

  data_ingrowth_ingrowth <- select(data_ingrowth, plotID, year, ingrowth_small, ingrowth_big) %>%
    group_by(plotID, year) %>% summarise_all(.funs = sum, na.rm = TRUE)

  data_ingrowth <- merge(data_ingrowth_stand, data_ingrowth_ingrowth, by = c("plotID", "year"))

  ######################################################################################

  initial_df <- data



  # 1 Calculate heights
  initial_df <- height_prediction(df_fit = data_height, df_predict = initial_df,
                                  species_n_threshold = species_n_threshold,
                                  height_model = height_model,
                                  BRNN_neurons = BRNN_neurons,
                                  height_pred_level = height_pred_level,
                                  eval_model_height = FALSE)$data_height_predictions

  initial_df <- crownHeight_prediction(df_fit = data_crownHeight,
                                       df_predict = initial_df,
                                       site_vars = site_vars,
                                       crownHeight_model = crownHeight_model,
                                       BRNN_neurons = BRNN_neurons,
                                       species_n_threshold = species_n_threshold,
                                       k = k,
                                       eval_model_crownHeight = FALSE)$predicted_crownHeight

  # 1 Mortality
  # 2 (mortality), 0 (normal), 1 (harvested), 3 (ingrowth small), 15 (ingrowth big)
  # Before you can apply mortality, you should have predicted height and crown height!

  # calculate volume
  if (volume_calculation == "form_factors"){

    initial_df$volume <- initial_df$height * initial_df$BA * initial_df$form
    initial_df$p_volume <- initial_df$p_height * initial_df$p_BA * initial_df$form

    # sum(data$p_volume > data$volume, na.rm = T)

  } else if (volume_calculation == "volume_functions"){

    initial_df$volume <- NA
    initial_df$p_volume <- NA

    if (is.null(data_volF_param)){

      stop("data_volF_param is not provided")

    }

    initial_df <- V_general(df = initial_df, data_volF_param = data_volF_param)

  } else if (volume_calculation == "tariffs"){

    if (is.null(data_tariffs)){

      stop("data_tariffs is not supplied!")
    }

    initial_df$volume <- NA
    initial_df$p_volume <- NA

    initial_df <- vol_tariffs(df = initial_df, data_tariffs = data_tariffs)

  } else if(volume_calculation == "slo_2p_volume_functions"){

    initial_df$volume <- NA
    initial_df$p_volume <- NA

    if (merchantable_whole_tree == "merchantable"){
      initial_df <- volume_merchantable(df = initial_df)
    } else if (merchantable_whole_tree == "whole_tree"){
      initial_df <- volume_whole_tree(df = initial_df)
    }

  } else{

    stop("Please define volume calculations: form_factors, volume_functions or tariffs")

  }

  list_results <- list()

  # Simulation starts and we save initial_df (without dead and harvested trees)
  list_results[[1]] <- initial_df



  # I create this empty object in case of no evaluation
  eval_mortality_output <- "the argument set_eval_mortality is set to FALSE"
  eval_ingrowth_output <- "the argument set_eval_ingrowth is set to FALSE"
  eval_BAI_output <- "the argument set_eval_BAI is set to FALSE"

  # If the length of mortality_sahre is 1, we replicate the value using the sim_steps
  if (length(mortality_share )== 1){

    mortality_share <- rep(mortality_share, sim_steps)

  }

  if (length(harvesting_sum)== 1){

    harvesting_sum <- rep(harvesting_sum, sim_steps)

  }

  if (length(mortality_share) < sim_steps && length(mortality_share) > 1){

    n_missing <- sim_steps - length(mortality_share)

    mortality_share <- c(mortality_share, rep(mortality_share[length(mortality_share)], n_missing))

    if (sim_mortality == TRUE){
      warning("The last value in mortality_share vector is used for undefined simulation years.")

    }
  }


  if (length(harvesting_sum) < sim_steps && length(harvesting_sum) > 1){

    n_missing <- sim_steps - length(harvesting_sum)

    harvesting_sum <- c(harvesting_sum, rep(harvesting_sum[length(harvesting_sum)], n_missing))

    if (sim_harvesting == TRUE){
      warning("The last value in harvesting_sum vector is used for undefined simulation years.")

    }
  }

  # This is only due to organization of the next for loop
  sim_steps <- sim_steps + 1

  sim = 2



  for (sim in 2:sim_steps){

    mortality_outputs <- predict_mortality(df_fit = data_mortality,
                                           df_predict = initial_df,
                                           df_climate = data_climate,
                                           site_vars = site_vars,
                                           sim_mortality = sim_mortality,
                                           mortality_model = mortality_model,
                                           nb_laplace = nb_laplace,
                                           rf_mtry = rf_mtry,
                                           mortality_share = mortality_share[sim-1],
                                           include_climate = include_climate,
                                           select_months_climate = select_months_climate,
                                           eval_model_mortality = set_eval_mortality,
                                           k = k, blocked_cv = blocked_cv,
                                           sim_step_years = sim_step_years)

    initial_df <- mortality_outputs$predicted_mortality

    if (set_eval_mortality == TRUE){

      eval_mortality_output <- mortality_outputs$eval_mortality

      set_eval_mortality <- FALSE

    }

    # Simulate harvesting
    if (sim_harvesting == TRUE){
      initial_df <- simulate_harvesting(df = initial_df,
                                        harvesting_sum = harvesting_sum[sim-1],
                                        forest_area_ha = forest_area_ha,
                                        harvesting_type = harvesting_type,


                                        harvest_sum_level = harvest_sum_level,
                                        plot_upscale_type = plot_upscale_type,
                                        plot_upscale_factor = plot_upscale_factor,


                                        final_cut_weight = final_cut_weight,
                                        thinning_small_weight = thinning_small_weight
                                        )
    }


    # 2 BAI
    BAI_outputs <- BAI_prediction(df_fit = data_BAI,
                                 df_predict = initial_df,
                                 site_vars = site_vars,
                                 rf_mtry = rf_mtry,
                                 species_n_threshold = species_n_threshold,
                                 include_climate = include_climate,
                                 eval_model_BAI = set_eval_BAI,
                                 k = k, blocked_cv = blocked_cv)

    # Na You might lose trees without BAI measurements! Be aware. Include warning!

    initial_df <- BAI_outputs$predicted_BAI

    if (set_eval_BAI == TRUE){

      eval_BAI_output <- BAI_outputs$eval_BAI

      set_eval_BAI <- FALSE

    }


    # 3 Calculate BAL
    initial_df <- calculate_BAL(initial_df)

    # 4 Calculate stand variables
    initial_df <- calculate_standVars(df = initial_df)

    # 5 Ingrowth
    ingrowth_outputs <- predict_ingrowth(df_fit = data_ingrowth, df_predict = initial_df,
                                   site_vars = site_vars, form_factors = form_factors,
                                   eval_model_ingrowth = set_eval_ingrowth,
                                   k = k, blocked_cv = blocked_cv, ingrowth_model = ingrowth_model
                                   )

    initial_df <- ingrowth_outputs$predicted_ingrowth

    if (set_eval_ingrowth == TRUE){

      eval_ingrowth_output <- ingrowth_outputs$eval_ingrowth

      set_eval_ingrowth <- FALSE

    }

    # We again calculate BAL and stand variables to update those variables by considering ingrowth
    initial_df <- calculate_BAL(initial_df)

    # 4.1 stand
    initial_df <- calculate_standVars(df = initial_df)


    # 5 Impute heights
    initial_df <- height_prediction(df_fit = data_height, df_predict = initial_df,
                                    species_n_threshold = species_n_threshold,
                                    height_model = height_model,
                                    BRNN_neurons = BRNN_neurons,
                                    height_pred_level = height_pred_level,
                                    eval_model_height = FALSE)$data_height_predictions

    # 6 Impute crown heights
    initial_df <- crownHeight_prediction(df_fit = data_crownHeight,
                                         df_predict = initial_df,
                                         site_vars = site_vars,
                                         crownHeight_model = crownHeight_model,
                                         BRNN_neurons = BRNN_neurons,
                                         species_n_threshold = species_n_threshold,
                                         k = k,
                                         eval_model_crownHeight = FALSE)$predicted_crownHeight



    # 7 Calculate Volume
    if (volume_calculation == "form_factors"){

      initial_df$volume <- initial_df$height * initial_df$BA * initial_df$form

      initial_df$p_volume <- initial_df$p_height * initial_df$p_BA * initial_df$form

    } else if (volume_calculation == "volume_functions"){

      initial_df$volume <- NA
      initial_df$p_volume <- NA

      if (is.null(data_volF_param)){

        stop("data_volF_param is not provided")

      }

      initial_df <- V_general(df = initial_df, data_volF_param)

    } else if (volume_calculation == "tariffs"){

      if (is.null(data_tariffs)){

        stop("data_tariffs is not supplied!")
      }

      initial_df$volume <- NA
      initial_df$p_volume <- NA

      initial_df <- vol_tariffs(df = initial_df, data_tariffs = data_tariffs)

    } else if(volume_calculation == "slo_2p_volume_functions"){

      initial_df$volume <- NA
      initial_df$p_volume <- NA

      if (merchantable_whole_tree == "merchantable"){
        initial_df <- volume_merchantable(df = initial_df)
      } else if (merchantable_whole_tree == "whole_tree"){
        initial_df <- volume_whole_tree(df = initial_df)
      }

    } else {

      stop("Please define volume calculations: form_factors, volume_functions or tariffs")

    }

    # 8 Save results

    list_results[[sim]] <- initial_df

    setTxtProgressBar(pb,sim)



  }

  final_ouputs <- list(

    sim_results = do.call(rbind, list_results),
    height_eval = h_predictions$data_height_eval,
    crownHeight_eval = Crown_h_predictions$eval_crownHeight,
    mortality_eval = eval_mortality_output,
    ingrowth_eval = eval_ingrowth_output,
    BAI_eval = eval_BAI_output

  )

  close(pb) # just to make print possible warning messages into a new row

  return(final_ouputs)

}

