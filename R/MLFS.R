#' MLFS
#'
#' Machine Learning Forest Simulator
#'
#' @param data_NFI data frame with individual tree variables
#' @param data_site data frame with site descriptors. This data is related to
#' data_NFI based on the 'plotID' column
#' @param data_tariffs optional, but mandatory if volume is calculated using the
#' one-parametric tariff functions. Data frame with plotID, species and V45. See
#' details.
#' @param data_climate data frame with climate data, covering the initial
#' calibration period and all the years which will be included in the simulation
#' @param thinning_weights_species data frame with thinning weights for each
#' species. The first column represents species code, each next column consists
#' of species-specific thinning weights applied in each simulation step
#' @param final_cut_weights_species data frame with final cut weights for each
#' species. The first column represents species code, each next column consists
#' of species-specific final cut weights applied in each simulation step
#' @param thinning_weights_plot data frame with harvesting weights related to plot
#' IDs, used for thinning
#' @param final_cut_weights_plot data frame with harvesting weights related
#' to plot IDs, used for final cut
#' @param sim_mortality logical, should mortality be simulated?
#' @param sim_ingrowth logical, should ingrowth be simulated?
#' @param sim_crownHeight logical, should crown heights be simulated? If TRUE,
#' a crownHeight column is expected in data_NFI
#' @param data_volF_param optional, data frame with species-specific volume
#' function parameters
#' @param form_factors optional, data frame with species-specific form factors
#' @param form_factors_level character, the level of specified form factors. It
#' can be 'species', 'plot' or 'species_plot'
#' @param uniform_form_factor numeric, uniform form factor to be used for all
#' species and plots. Only if form_factors are not provided
#' @param volume_calculation character string defining the method for volume
#' calculation: 'tariffs', 'volume_functions', 'form_factors' or
#' 'slo_2p_volume_functions'
#' @param sim_harvesting logical, should harvesting be simulated?
#' @param harvest_sum_level integer with value 0 or 1 defining the level of
#' specified harvesting sum: 0 for plot level and 1 for regional level.
#' @param plot_upscale_type character defining the upscale method of plot level
#' area. It can be 'area' or 'upscale factor'. If 'area', provide the forest
#' area represented by all plots in hectares (forest_area_ha argument). If
#' 'factor', provide the fixed factor to upscale the area of all plots. Please
#' note: forest_area_ha/plot_upscale_factor = number of unique plots. This
#' argument is important when harvesting sum is defined on regional level.
#' @param plot_upscale_factor numeric value to be used to upscale area of each
#' plot
#' @param harvesting_sum a value, or a vector of values defining the harvesting
#' sums through the simulation stage. If a single value, then it is used in all
#' simulation steps. If a vector of values, the first value is used in the first
#' step, the second in the second step, etc.
#' @param sim_mortality logical, should mortality be simulated?
#' @param mortality_share a value, or a vector of values defining the proportion
#' of the volume which is to be the subject of mortality. If a single value,
#' then it is used in all simulation steps. If a vector of values, the first
#' value is used in the first step, the second in the second step, and so on.
#' @param mortality_share_type character, it can be 'volume' or 'n_trees'. If
#' 'volume' then the mortality share relates to total standing volume, if
#' 'n_trees' then mortality share relates to the total number of standing trees
#' @param forest_area_ha the total area of all forest which are subject of the
#' simulation.
#' @param harvesting_type character, it could be 'random', 'final_cut',
#' 'thinning' or 'combined'. The latter combines 'final_cut' and 'thinning'
#' options, where the share of each is specified with the argument
#' 'share_thinning'
#' @param share_thinning numeric, a number oa a vector of numbers between 0 and
#' 1 that specifies the share of thinning in comparison to final_cut. Only used
#' if harvesting_type is 'combined'
#' @param final_cut_weight numeric value affecting the probability distribution
#' of harvested trees. Greater value increases the share of harvested trees
#' having larger DBH. Default is 10.
#' @param thinning_small_weight numeric value affecting the probability
#' distribution of harvested trees. Greater value increases the share of
#' harvested trees having smaller DBH. Default is 1.
#' @param k the number of folds to be used in the k fold cross-validation
#' @param blocked_cv logical, should the blocked cross-validation be used in the
#' evaluation phase?
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
#' generalized linear models (Poisson regression), 'ZIF_poiss' for zero inflated
#' Poisson regression and 'rf' for random forest
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
#' @param max_size a data frame with the maximum values of DBH for each species.
#' If a tree exceeds this value, it dies. If not provided, the maximum is
#' estimated from the input data. Two columns must be present, i.e. 'species'
#' and 'DBH_max'
#' @param max_size_increase_factor numeric value, which will be used to increase
#' the max DBH for each species, when the maximum is estimated from the input
#' data. If the argument 'max_size' is provided, the 'max_size_increase_factor'
#' is ignored. Default is 1. To increase maximum for 10 percent, use 1.1.
#' @param ingrowth_codes numeric value or a vector of codes which refer to
#' ingrowth trees
#' @param ingrowth_max_DBH_percentile which percentile should be used to estimate
#' the maximum simulated value of ingrowth trees?
#' @param measurement_thresholds data frame with two variables: 1) DBH_threshold
#' and 2) weight. This information is used to assign the correct weights in BAI
#' and increment sub-model; and to upscale plot-level data to hectares.
#' @param area_correction optional data frame with three variables: 1) plotID and
#' 2) DBH_threshold and 3) the correction factor to be multiplied by weight for
#' this particular category.
#' @param export_csv logical, if TRUE, at each simulation step, the results are
#' saved in the current working directory as csv file
#' @param sim_export_mode logical, if FALSE, the results of the individual
#' simulation steps are not merged into the final export table. Therefore,
#' output element 1 ($sim_results) will be empty. This was introduced to allow
#' simulations when using larger data sets and long term simulations that might
#' exceed the available RAM. In such cases, we recommend setting the argument
#' export_csv = TRUE, which will export each simulation step to the current
#' working directory.
#' @param include_mortality_BAI logical, should basal area increments (BAI) be
#' used as independent variable for predicting individual tree morality?
#' @param intermediate_print logical, if TRUE intermediate steps will be printed
#' while MLFS is running
#'
#' @examples
#' \dontrun{
#' library(MLFS)
#'
#' # open example data
#' data(data_NFI)
#' data(data_site)
#' data(data_climate)
#' data(data_volF_param)
#' data(measurement_thresholds)
#'
#' test_simulation <- MLFS(data_NFI = data_NFI,
#'  data_site = data_site,
#'  data_climate = data_climate,
#'  data_volF_param = data_volF_param,
#'  form_factors = form_factors,
#'  sim_steps = 2,
#'  sim_harvesting = TRUE,
#'  harvesting_sum = 100000,
#'  harvest_sum_level = 1,
#'  plot_upscale_type = "factor",
#'  plot_upscale_factor = 1600,
#'  measurement_thresholds = measurement_thresholds,
#'  ingrowth_codes = c(3,15),
#'  volume_calculation = "volume_functions",
#'  select_months_climate = seq(6,8),
#'  intermediate_print = FALSE
#'  )
#' }

MLFS <- function(data_NFI, data_site,
                 data_tariffs = NULL,
                 data_climate = NULL,
                 data_volF_param = NULL,

                 thinning_weights_species = NULL,
                 final_cut_weights_species = NULL,
                 thinning_weights_plot = NULL,
                 final_cut_weights_plot = NULL,

                 form_factors = NULL,
                 form_factors_level = 'species_plot',
                 uniform_form_factor = 0.42,
                 sim_steps,
                 volume_calculation = "volume_functions",
                 merchantable_whole_tree = "merchantable",

                 sim_harvesting = TRUE,
                 sim_mortality = TRUE,
                 sim_ingrowth = TRUE,
                 sim_crownHeight = TRUE,

                 harvesting_sum = NULL,
                 forest_area_ha = NULL,
                 harvest_sum_level = NULL,
                 plot_upscale_type = NULL,
                 plot_upscale_factor = NULL,

                 mortality_share = NA,
                 mortality_share_type = "volume",
                 mortality_model = "glm",
                 ingrowth_model = "ZIF_poiss",
                 rf_mtry = NULL,
                 nb_laplace = 0,
                 harvesting_type = "final_cut",
                 share_thinning = 0.80,
                 final_cut_weight = 10,
                 thinning_small_weight = 1,

                 species_n_threshold = 100,
                 height_model = "brnn",
                 crownHeight_model = "brnn",
                 BRNN_neurons = 3,
                 height_pred_level = 0, # prediction level for lmfor (0 regional, 1 plot level)
                 include_climate = FALSE,
                 select_months_climate = c(1,12),
                 set_eval_mortality = TRUE,
                 set_eval_crownHeight = TRUE,
                 set_eval_height = TRUE,
                 set_eval_ingrowth = TRUE,
                 set_eval_BAI = TRUE,
                 k = 10, blocked_cv = TRUE,
                 max_size = NULL,
                 max_size_increase_factor = 1.0,
                 ingrowth_codes = c(3),
                 ingrowth_max_DBH_percentile = 0.90,
                 measurement_thresholds = NULL,
                 area_correction = NULL,
                 export_csv = FALSE,
                 sim_export_mode = TRUE,
                 include_mortality_BAI = TRUE,
                 intermediate_print = FALSE
                 ){

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
  n <- NULL
  species <- NULL
  DBH_max <- NULL
  max_DBH_data <- NULL
  max_size_DBH_joint <- NULL
  DBH_max_data <- NULL
  weight <- NULL
  area_factor <- NULL

  # NFI codes
  ## 0  (normal)
  ## 1  (harvested)
  ## 2  (dead - mortality)
  ## 3  (ingrowth level 1) # optional
  ## 15 (ingrowth level 2) # optional

  pb = txtProgressBar(min = 0, max = sim_steps, initial = 0, style = 3)

  setTxtProgressBar(pb, 1)

  options(dplyr.summarise.inform= FALSE)

  # calculate maximum tree size
  max_size_data <- dplyr::group_by(data_NFI, species) %>% summarise(DBH_max_data = max(DBH, na.rm = TRUE))

  # Increase
  max_size_data$DBH_max_data <- max_size_data$DBH_max_data * max_size_increase_factor

  # The measurement_threshold table must be specified
  if (is.null(measurement_thresholds)){

    stop(paste0("measurement_thresholds table is missing. This is a data frame ",
                "with two variables: 1) DBH_threshold and 2) weight. This " ,
                "information is used to assign the correct weights in BAI and ",
                "increment sub-model; and to upscale plot-level data to hectares."))

  }

  # If not provided by user, we use the calculations
  if (!is.null(max_size)){

    if (sum(colnames(max_size) %in% c('species', "DBH_max")) < 2){

      stop(paste0("max_DBH data frame should have two columns 'species' and 'DBH_max'"))
    }

    max_size_data <- merge(max_size_data, max_size, by = "species", all.x = TRUE)
    max_size_data <-  mutate(max_size_data,
                             max_size_DBH_joint = ifelse(is.na(DBH_max), DBH_max_data, DBH_max),
                             DBH_max = NULL,
                             DBH_max_data= NULL) %>%
      mutate(BA_max = ((max_size_DBH_joint/2)^2 * pi)/10000,
             max_size_DBH_joint = NULL
             )

  } else {

    max_size_data <- mutate(max_size_data,
                       BA_max = ((DBH_max_data/2)^2 * pi)/10000,
                       DBH_max_data = NULL)
  }

  # String related to height model is converted to lowercase
  height_model <- tolower(height_model)
  crownHeight_model <- tolower(crownHeight_model)

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
  site_vars <- colnames(data_site)[!(colnames(data_site) %in% c("plotID"))]

  # merge NFI and site descriptors
  data <- merge(data_NFI, data_site, by = "plotID")

  # remove the objects to make free space
  rm(data_NFI)
  rm(data_site)

  # Function to calculate the most common value in a vector
  # source: https://stackoverflow.com/questions/29255473/most-frequent-value-mode-by-group
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  # 1 Calculate basal area and remove DBH, we don't need it anymore
  data <- dplyr::mutate(data, BA = ((DBH/2)^2 * pi)/10000)

  # Calculate measurement thresholds in terms of basal area
  measurement_thresholds$BA_threshold <- ((measurement_thresholds$DBH_threshold/2)^2 * pi)/10000

  data <- data %>% mutate(weight = ifelse(BA >= max(measurement_thresholds$BA_threshold),
                                          measurement_thresholds[, "weight"][which.max(measurement_thresholds$BA_threshold)],
                                          measurement_thresholds[, "weight"][which.min(measurement_thresholds$BA_threshold)]),

                          DBH_threshold = ifelse(BA >= max(measurement_thresholds$BA_threshold),
                                                 measurement_thresholds[, "DBH_threshold"][which.max(measurement_thresholds$BA_threshold)],
                                                 measurement_thresholds[, "DBH_threshold"][which.min(measurement_thresholds$BA_threshold)]))

  # In case area correction factors are provided we use them to correct plot weights
  if (!is.null(area_correction)){

    data <- merge(data, area_correction, by = c("plotID", "DBH_threshold"), all.x = TRUE)

    data <- dplyr::mutate(data, area_factor = ifelse(is.na(area_factor), 1, area_factor),
                          weight = weight*area_factor, area_factor = NULL, DBH_threshold = NULL)
  } else {

    area_correction <- NULL
    data$DBH_threshold <- NULL

  }

  if (sim_ingrowth == TRUE){

  # calculate the ingrowth table
  ingrowth_table <- dplyr::filter(data, code %in% c(ingrowth_codes)) %>%
    dplyr::group_by(code) %>%
    dplyr::summarise(DBH_threshold = min(DBH),
              DBH_max = quantile(DBH, ingrowth_max_DBH_percentile, na.rm = T),
              weight = Mode(weight))

  # calculate the parameters for ingrowth distributions
  ing_param_list <- list()
  ipl_holder <- 1

  for (i in unique(ingrowth_table$code)){

    temp_DBH_max <- dplyr::filter(ingrowth_table, code == i) %>%
      dplyr::select(DBH_max)

    temp_par_table <- dplyr::filter(data, code == i,
                                    DBH < as.numeric(temp_DBH_max))

    temp_parameters <- quantile(temp_par_table$DBH, probs = seq(0, 1, 0.05))

    ing_param_list[[ipl_holder]] <- temp_parameters
    names(ing_param_list)[ipl_holder] <- i
    ipl_holder <- ipl_holder + 1

  }
}

  data$DBH <- NULL
  data <- add_stand_variables(df = data)

  # 2 Calculate BAL
  data$BAL <- NA
  data <- calculate_BAL(df = data)


  if (sim_crownHeight == TRUE){

    # If crown heights will be simulated, the crownHeight column is expected in data_NFI
    if (!("crownHeight" %in% colnames(data))){

      stop("crownHeight variable must be present in data_NFI")

    }

  } else {

    data$crownHeight <- NA

  }

  # 3 Transform data, where yield is properly expressed
  data <- transform_data(df = data, include_climate = include_climate,
                         df_climate = data_climate, select_months_climate = select_months_climate
                         )

  ######################################################
  # create fitting data frames for
  # 1) height
  data_height <- dplyr::filter(data, !is.na(height))

  # 2) crown height
  if (sim_crownHeight == TRUE){
    data_crownHeight <- dplyr::filter(data, !is.na(crownHeight))
  }

  # 3) BAI
  data_BAI <- dplyr::filter(data, !is.na(BAI))

  # 4) mortality
  data_mortality <- data

  # This part needs to be run to get the elevation results
  h_predictions <- height_prediction(df_fit = data_mortality,
                                   df_predict = data_mortality,
                                   species_n_threshold = species_n_threshold,
                                   height_pred_level = height_pred_level,
                                   height_model = height_model,
                                   BRNN_neurons = BRNN_neurons,
                                   eval_model_height = set_eval_height,
                                   blocked_cv = blocked_cv, k = k)

  data_mortality <- h_predictions$data_height_predictions

  if (sim_crownHeight == TRUE & sum(is.na(data_mortality$crownHeight)) > 0){

    Crown_h_predictions <- crownHeight_prediction(df_fit = data_mortality,
                                                  df_predict = data_mortality,
                                                  site_vars = site_vars,
                                                  crownHeight_model = crownHeight_model,
                                                  BRNN_neurons = BRNN_neurons,
                                                  species_n_threshold = species_n_threshold,
                                                  k = k, blocked_cv = blocked_cv,
                                                  eval_model_crownHeight = set_eval_crownHeight)

    data_mortality <- Crown_h_predictions$predicted_crownHeight

  } else {

    Crown_h_predictions <- list()
    Crown_h_predictions$eval_crownHeight <- "crownHeight is not simulated"
    Crown_h_predictions$model_species <- "crownHeight is not simulated"
    Crown_h_predictions$model_speciesGroups <- "crownHeight is not simulated"

  }

  data_mortality <- dplyr::filter(data_mortality, !is.na(BA))

  # 5) Ingrowth - consists only of tree codes 0, 3, 15
   if (sim_ingrowth == TRUE){

     data_ingrowth <- data_mortality

     ing_codes <- unique(ingrowth_table$code)
     ing_code_var <- c() # This is empty vector where we store ingrowth_var names

     data_ingrowth <- dplyr::filter(data_ingrowth, code %in% c(0, ing_codes))

     # Here we define different levels of ingrowth
     for (i_codes in ing_codes){

       var_name <- paste0("ingrowth_", i_codes)
       ing_code_var <- c(ing_code_var, var_name)

       data_ingrowth$new_var <- ifelse(data_ingrowth$code == i_codes, 1, 0)
       colnames(data_ingrowth)[ncol(data_ingrowth)] <- var_name

     }

     data_ingrowth_stand <- dplyr::select(data_ingrowth, plotID, year, stand_BA, stand_n, BAL, all_of(site_vars)) %>%
       group_by(plotID, year) %>% summarise_all(.funs = mean, na.rm = TRUE)

     data_ingrowth_ingrowth <- dplyr::select(data_ingrowth, plotID, year, all_of(ing_code_var)) %>%
       group_by(plotID, year) %>% summarise_all(.funs = sum, na.rm = TRUE)

     data_ingrowth <- merge(data_ingrowth_stand, data_ingrowth_ingrowth, by = c("plotID", "year"))

   }

  ######################################################################################

  initial_df <- data

  # remove data
  rm(data)

  # 1 Calculate heights
  if (sum(is.na(initial_df$height)) > 0){

  initial_df <- height_prediction(df_fit = data_height, df_predict = initial_df,
                                  species_n_threshold = species_n_threshold,
                                  height_model = height_model,
                                  BRNN_neurons = BRNN_neurons,
                                  height_pred_level = height_pred_level,
                                  eval_model_height = FALSE)$data_height_predictions

  }


  # 2 Calculate CrownHeights
  if (sim_crownHeight == TRUE & sum(is.na(initial_df$crownHeight)) > 0){

  initial_df <- crownHeight_prediction(df_fit = data_crownHeight,
                                       df_predict = initial_df,
                                       site_vars = site_vars,
                                       crownHeight_model = crownHeight_model,
                                       BRNN_neurons = BRNN_neurons,
                                       species_n_threshold = species_n_threshold,
                                       k = k,
                                       eval_model_crownHeight = FALSE)$predicted_crownHeight
  }

  # 3 calculate volume
  if (volume_calculation == "form_factors"){

    initial_df$volume <- NA
    initial_df$p_volume <- NA

    initial_df <- vol_form_factors(df = initial_df, form_factors = form_factors,
                                   form_factors_level = form_factors_level,
                                   uniform_form_factor = uniform_form_factor)

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


  # For all trees, simulated BAI for half of the period
  # This is crucial in terms of correct harvesting and mortality estimates

  # Simulate BAI for halfPeriod

  # sort(table(initial_df$species))

  initial_df <- BAI_prediction_halfPeriod(df_fit = data_BAI,
                                          df_predict = initial_df,
                                          site_vars = site_vars,
                                          rf_mtry = rf_mtry,
                                          species_n_threshold = species_n_threshold,
                                          include_climate = include_climate,
                                          measurement_thresholds = measurement_thresholds,
                                          area_correction = area_correction
                                          )

  # Next, we simulate height and crownHeight based on half period attributes
  # Calculate tree heights - half Period

  initial_df <- height_prediction_halfPeriod(df_fit = data_height, df_predict = initial_df,
                                             species_n_threshold = species_n_threshold,
                                             height_model = height_model,
                                             BRNN_neurons = BRNN_neurons,
                                             height_pred_level = height_pred_level)

  # Calculate tree crownHeights half Period

 if (sim_crownHeight == TRUE ){

  initial_df <- crownHeight_prediction_halfPeriod(df_fit = data_crownHeight,
                                                  df_predict = initial_df,
                                                  site_vars = site_vars,
                                                  crownHeight_model = crownHeight_model,
                                                  BRNN_neurons = BRNN_neurons,
                                                  species_n_threshold = species_n_threshold)
 } else {

   initial_df$crownHeight_mid <- NA

 }

  # Calculate Volume
  if (volume_calculation == "form_factors"){

    initial_df <- vol_form_factors_halfPeriod(df = initial_df, form_factors = form_factors,
                                              form_factors_level = form_factors_level,
                                              uniform_form_factor = uniform_form_factor)

  } else if (volume_calculation == "volume_functions"){

    initial_df$volume_mid <- NA

    if (is.null(data_volF_param)){

      stop("data_volF_param is not provided")

    }

    initial_df <- V_general_halfPeriod(df = initial_df, data_volF_param)

  } else if (volume_calculation == "tariffs"){

    if (is.null(data_tariffs)){

      stop("data_tariffs is not supplied!")
    }

    initial_df$volume_mid <- NA

    initial_df <- vol_tariffs_halfPeriod(df = initial_df, data_tariffs = data_tariffs)

  } else if(volume_calculation == "slo_2p_volume_functions"){

    initial_df$volume_mid <- NA

    if (merchantable_whole_tree == "merchantable"){
      initial_df <- volume_merchantable_halfPeriod(df = initial_df)
    } else if (merchantable_whole_tree == "whole_tree"){
      initial_df <- volume_whole_tree_halfPeriod(df = initial_df)
    }

  } else {

    stop("Please define volume calculations: form_factors, volume_functions or tariffs")

  }

  # Ingrowth is now added, BAL and stand variables have changed. We therefore
  # update all variables
  initial_df <- calculate_BAL_halfPeriod(df = initial_df)
  initial_df <- add_stand_variables_halfPeriod(df = initial_df)

  list_results <- list()

  initial_df$p_BA_mid <- NA
  initial_df$p_weight_mid <- NA
  initial_df$p_height_mid <- NA
  initial_df$p_crownHeight_mid <- NA
  initial_df$p_volume_mid <- NA

  # Simulation starts and we save initial_df (without dead and harvested trees)

  if (sim_export_mode == TRUE){

    list_results[[1]] <- initial_df

  }

  if (export_csv == TRUE){
    write.csv(initial_df, paste0("output_step_0.csv"), row.names = FALSE)
  }


  # I create this empty objects in case of no evaluation
  eval_mortality_output <- "the argument set_eval_mortality is set to FALSE"
  eval_ingrowth_output <- "the argument set_eval_ingrowth is set to FALSE"
  eval_BAI_output <- "the argument set_eval_BAI is set to FALSE"

  # If ingrowth is simulated, then these objects are later overwritten

  if (sim_ingrowth == TRUE){

    for (i_codes in ing_codes){
      assign(paste0("ing_model_output_", i_codes), "Ingrowth is not simulated. No model output available")
    }

  }

  # If the length of share_thinning is 1, we replicate the value using the sim_steps
  if (length(share_thinning) == 1){
    share_thinning <- rep(share_thinning, sim_steps)
  }

  # If the length of mortality_share is 1, we replicate the value using the sim_steps
  if (length(mortality_share) == 1){
    mortality_share <- rep(mortality_share, sim_steps)
  }

  # If the length of thinning_small_weight is 1, we replicate the value using the sim_steps
  if (length(thinning_small_weight) == 1){
    thinning_small_weight <- rep(thinning_small_weight, sim_steps)
  }

  # If the length of final_cut_weight is 1, we replicate the value using the sim_steps
  if (length(final_cut_weight) == 1){
    final_cut_weight <- rep(final_cut_weight, sim_steps)
  }

  # If the length of harvesting_sum is 1, we replicate the value using the sim_steps
  if (length(harvesting_sum)== 1){
    harvesting_sum <- rep(harvesting_sum, sim_steps)
  }

  # If the ncol of thinning_weights_species is 2, we replicate the column using the sim_steps
  if (!is.null(thinning_weights_species)){
    if (ncol(thinning_weights_species) == 2){

      thinning_weights_species_column <- thinning_weights_species[,2]

      for (missing_step in 2:sim_steps){

        thinning_weights_species[,missing_step +1] <- thinning_weights_species_column

      }
    }
  }

  # If the ncol of final_cut_weights_species is 2, we replicate the column using the sim_steps
  if (!is.null(final_cut_weights_species)){
    if (ncol(final_cut_weights_species) == 2){

      final_cut_weights_species_column <- final_cut_weights_species[,2]

      for (missing_step in 2:sim_steps){

        final_cut_weights_species[,missing_step +1] <- final_cut_weights_species_column

      }
    }
  }

  # If the ncol of final_cut_weights_plot is 2, we replicate the column using the sim_steps
  if (!is.null(final_cut_weights_plot)){
    if (ncol(final_cut_weights_plot) == 2){

      thinning_weights_plot_column <- final_cut_weights_plot[,2]

      for (missing_step in 2:sim_steps){

        final_cut_weights_plot[,missing_step +1] <- thinning_weights_plot_column

      }
    }
  }

  # If the ncol of thinning_weights_plot is 2, we replicate the column using the sim_steps
  if (!is.null(thinning_weights_plot)){
    if (ncol(thinning_weights_plot) == 2){

      thinning_weights_plot_column <- thinning_weights_plot[,2]

      for (missing_step in 2:sim_steps){

        thinning_weights_plot[,missing_step +1] <- thinning_weights_plot_column

      }
    }
  }

  if (length(mortality_share) < sim_steps && length(mortality_share) > 1){

    n_missing <- sim_steps - length(mortality_share)
    mortality_share <- c(mortality_share, rep(mortality_share[length(mortality_share)], n_missing))

    if (sim_mortality == TRUE){
      warning("The last value in mortality_share vector is used for undefined simulation years.")

    }
  }

  if (length(thinning_small_weight) < sim_steps && length(thinning_small_weight) > 1){

    n_missing <- sim_steps - length(thinning_small_weight)
    thinning_small_weight <- c(thinning_small_weight, rep(thinning_small_weight[length(thinning_small_weight)], n_missing))

    if (sim_mortality == TRUE){
      warning("The last value in thinning_small_weight vector is used for undefined simulation years.")

    }
  }

  if (length(final_cut_weight) < sim_steps && length(final_cut_weight) > 1){

    n_missing <- sim_steps - length(final_cut_weight)
    final_cut_weight <- c(final_cut_weight, rep(final_cut_weight[length(final_cut_weight)], n_missing))

    if (sim_mortality == TRUE){
      warning("The last value in final_cut_weight vector is used for undefined simulation years.")

    }
  }

  if (length(harvesting_sum) < sim_steps && length(harvesting_sum) > 1){

    n_missing <- sim_steps - length(harvesting_sum)
    harvesting_sum <- c(harvesting_sum, rep(harvesting_sum[length(harvesting_sum)], n_missing))

    if (sim_harvesting == TRUE){
      warning("The last value in harvesting_sum vector is used for undefined simulation years.")
    }
  }

  if (length(share_thinning) < sim_steps && length(share_thinning) > 1){

    n_missing <- sim_steps - length(share_thinning)
    share_thinning <- c(share_thinning, rep(share_thinning[length(share_thinning)], n_missing))

    if (sim_harvesting == TRUE){
      warning("The last value in share_thinning vector is used for undefined simulation years.")
    }
  }

  # If the ncol of thinning_weights_species is > 2, we replicate the column using the sim_steps
  if (!is.null(thinning_weights_species)){
    if (ncol(thinning_weights_species) > 2 && ncol(thinning_weights_species) < (sim_steps + 1)){

      thinning_weights_species_column <- thinning_weights_species[,2]

      for (missing_step in (ncol(thinning_weights_species):sim_steps)){

        thinning_weights_species[,missing_step +1] <- thinning_weights_species_column

      }
    }
  }

  # If the ncol of final_cut_weights_species is > 2, we replicate the column using the sim_steps
  if (!is.null(final_cut_weights_species)){
    if (ncol(final_cut_weights_species) > 2 && ncol(final_cut_weights_species) < (sim_steps + 1)){

      final_cut_weights_species_column <- final_cut_weights_species[,2]

      for (missing_step in (ncol(final_cut_weights_species):sim_steps)){

        final_cut_weights_species[,missing_step +1] <- final_cut_weights_species_column

      }
    }
  }

  # If the ncol of final_cut_weights_plot is > 2, we replicate the column using the sim_steps
  if (!is.null(final_cut_weights_plot)){
    if (ncol(final_cut_weights_plot) > 2 && ncol(final_cut_weights_plot) < (sim_steps + 1)){

      thinning_weights_plot_column <- final_cut_weights_plot[,2]

      for (missing_step in (ncol(final_cut_weights_plot):sim_steps)){

        final_cut_weights_plot[,missing_step +1] <- thinning_weights_plot_column

      }
    }
  }

  # If the ncol of thinning_weights_plot is > 2, we replicate the column using the sim_steps
  if (!is.null(thinning_weights_plot)){

    if (ncol(thinning_weights_plot) > 2 && ncol(thinning_weights_plot) < (sim_steps + 1)){

      thinning_weights_plot_column <- thinning_weights_plot[,2]

      for (missing_step in (ncol(thinning_weights_plot):sim_steps)){

        thinning_weights_plot[,missing_step +1] <- thinning_weights_plot_column

      }
    }
  }

  # If sim_harvesting = TRUE, harvesting_sum, harvest_sum_level, plot_upscale_tpye
  # and plot_upscale factor must be defined

  if (sim_harvesting == TRUE){

    if (is.null(harvesting_sum)){

      stop(paste0("The sim_harvesting is set to TRUE, but the harvesting volume is not specified.",
                 " Define the harveting volume with the harvesting_sum argument."))

    }


    if (is.null(harvesting_sum)){

      stop(paste0("The argument sim_harvesting is set to TRUE, but the harvesting_sum is not specified.",
                  "Please define the amount of harvested volume for each simulation step."))
      }

    if (is.null(harvest_sum_level)){

      stop(paste0("The argument sim_harvesting is set to TRUE, and the harvesting_sum is specified. ",
           "However, you should also specify the argument harvest_sum_level. ",
           "If the amount of harvested volume is defined on a plot level, set harvest_sum_level = 0, ",
           " If the amount of harvested volume is defined on a regional level, i.e. for all plots using the ",
           "upscale factor, set harvest_sum_level = 1."))

    }

    if ((harvest_sum_level == 1) & is.null(plot_upscale_type)){

      stop(paste0("The harvest_sum_level = 1 (harvesting is defined on regional level), ",
                  "but the plot_upscale_type is not defined. Please define the plot_upscale_type. ",
                  "It can be 'area' or 'upscale factor'." ))
    }


    if ((plot_upscale_type == 'area') &  is.null(forest_area_ha)){

      stop(paste0("plot_upscale_type is set to 'area', but the forest_area_ha is not defined. ",
                  "forest_area_ha is the total area of all forest which are subject of a simulation."))
    }

    if ((plot_upscale_type == 'factor') &  is.null(plot_upscale_factor)){

      stop(paste0("plot_upscale_type is set to 'factor', but the plot_upscale_factor is not defined. ",
                  "plot_upscale_factor is a value to be used to upscale area from plot to regional level."))
    }

  }

  # This is only due to organization of the next for loop
  sim_steps <- sim_steps + 1

  for (sim in 2:sim_steps){

    if (intermediate_print == TRUE){

      print(paste0("simulating mortality in step ", sim - 1))

    }

    # Simulate mortality
    mortality_outputs <- predict_mortality(df_fit = data_mortality,
                                           df_predict = initial_df,
                                           mortality_share_type = mortality_share_type,
                                           df_climate = data_climate,
                                           site_vars = site_vars,
                                           sim_mortality = sim_mortality,
                                           mortality_model = mortality_model,
                                           nb_laplace = nb_laplace,
                                           rf_mtry = rf_mtry, sim_crownHeight = sim_crownHeight,
                                           mortality_share = mortality_share[sim-1],
                                           include_climate = include_climate,
                                           select_months_climate = select_months_climate,
                                           eval_model_mortality = set_eval_mortality,
                                           k = k, blocked_cv = blocked_cv,
                                           sim_step_years = sim_step_years,
                                           df_max_size = max_size_data,
                                           ingrowth_codes = ingrowth_codes,
                                           include_mortality_BAI = include_mortality_BAI,
                                           intermediate_print = intermediate_print
                                           )

    initial_df <- mortality_outputs$predicted_mortality

    # Save the model for the final list
    mortality_output_model <- mortality_outputs$model_output

    if (set_eval_mortality == TRUE){

      eval_mortality_output <- mortality_outputs$eval_mortality

      set_eval_mortality <- FALSE

    }

    #remove mortality outputs to clear space
    rm(mortality_outputs)

    # Simulate harvesting
    if (sim_harvesting == TRUE){

      if (intermediate_print == TRUE){

        print(paste0("simulating harvesting in step ", sim - 1))

      }



      initial_df <- simulate_harvesting(df = initial_df,
                                        harvesting_sum = harvesting_sum[sim-1],
                                        forest_area_ha = forest_area_ha,
                                        harvesting_type = harvesting_type,
                                        share_thinning = share_thinning[sim-1],

                                        df_thinning_weights_species = if (!is.null(thinning_weights_species)) df_weights <- thinning_weights_species[,c(1,sim)],
                                        df_final_cut_weights_species = if (!is.null(final_cut_weights_species)) df_weights <- final_cut_weights_species[,c(1,sim)],

                                        df_thinning_weights_plot = if (!is.null(thinning_weights_plot)) df_weights <- thinning_weights_plot[,c(1,sim)],
                                        df_final_cut_weights_plot = if (!is.null(final_cut_weights_plot)) df_weights <- final_cut_weights_plot[,c(1,sim)],

                                        harvest_sum_level = harvest_sum_level,
                                        plot_upscale_type = plot_upscale_type,
                                        plot_upscale_factor = plot_upscale_factor,

                                        final_cut_weight = final_cut_weight[sim-1],
                                        thinning_small_weight = thinning_small_weight[sim-1])
    }

    # Simulate BAI

    if (intermediate_print == TRUE){

      print(paste0("simulating BAI in step ", sim - 1))

    }

    BAI_outputs <- BAI_prediction(df_fit = data_BAI,
                                 df_predict = initial_df,
                                 site_vars = site_vars,
                                 rf_mtry = rf_mtry,
                                 species_n_threshold = species_n_threshold,
                                 include_climate = include_climate,
                                 eval_model_BAI = set_eval_BAI,
                                 k = k, blocked_cv = blocked_cv,
                                 measurement_thresholds = measurement_thresholds,
                                 area_correction = area_correction)

    BAI_outputs_model_species <- BAI_outputs$rf_model_species
    BAI_outputs_model_groups <- BAI_outputs$rf_model_speciesGroups

    # You might lose trees without BAI measurements! Be aware. Include warning!
    initial_df <- BAI_outputs$predicted_BAI

    if (set_eval_BAI == TRUE){

      eval_BAI_output <- BAI_outputs$eval_BAI

      set_eval_BAI <- FALSE

    }

    # remove the object to save space
    remove(BAI_outputs)

    # Calculate BAL
    initial_df <- calculate_BAL(initial_df)

    # Calculate stand variables
    initial_df <- add_stand_variables(df = initial_df)

    # Simulate Ingrowth
    if (sim_ingrowth == TRUE){

      if (intermediate_print == TRUE){

        print(paste0("simulating ingrowth in step ", sim - 1))

      }

    ingrowth_outputs <- predict_ingrowth(df_fit = data_ingrowth, df_predict = initial_df,
                                   site_vars = site_vars, include_climate = include_climate,
                                   eval_model_ingrowth = set_eval_ingrowth,
                                   k = k, blocked_cv = blocked_cv, ingrowth_model = ingrowth_model,
                                   ingrowth_table = ingrowth_table,
                                   DBH_distribution_parameters = ing_param_list)

    initial_df <- ingrowth_outputs$predicted_ingrowth

    # save models for ingrowth
    for (i_code in ing_codes){

      assign(paste0("ing_model_output_", i_code), ingrowth_outputs[[paste0("mod_ingrowth_", i_code)]])

    }

    if (set_eval_ingrowth == TRUE){

      eval_ingrowth_output <- ingrowth_outputs$eval_ingrowth

      set_eval_ingrowth <- FALSE

    }

    # remove the ingrowth_outputs to save space
    rm(ingrowth_outputs)

    # Ingrowth is now added, BAL and stand variables have changed. We therefore
    # update all variables
    initial_df <- calculate_BAL(initial_df)
    initial_df <- add_stand_variables(df = initial_df)

    }

    # Calculate tree heights

    if (intermediate_print == TRUE){

      print(paste0("updating tree heights in step ", sim - 1))

    }

    initial_df <- height_prediction(df_fit = data_height, df_predict = initial_df,
                                    species_n_threshold = species_n_threshold,
                                    height_model = height_model,
                                    BRNN_neurons = BRNN_neurons,
                                    height_pred_level = height_pred_level,
                                    eval_model_height = FALSE)$data_height_predictions

    # Calculate tree crownHeights
    if (sim_crownHeight == TRUE){

      if (intermediate_print == TRUE){

        print(paste0("updating crown heights in step ", sim - 1))

      }

    initial_df <- crownHeight_prediction(df_fit = data_crownHeight,
                                         df_predict = initial_df,
                                         site_vars = site_vars,
                                         crownHeight_model = crownHeight_model,
                                         BRNN_neurons = BRNN_neurons,
                                         species_n_threshold = species_n_threshold,
                                         k = k,
                                         eval_model_crownHeight = FALSE)$predicted_crownHeight
    }


    if (intermediate_print == TRUE){

      print(paste0("calculating tree volume in step ", sim - 1))

    }

    # Calculate Volume
    if (volume_calculation == "form_factors"){

      initial_df <- vol_form_factors(df = initial_df, form_factors = form_factors,
                                     form_factors_level = form_factors_level,
                                     uniform_form_factor = uniform_form_factor)

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

    # For all trees, simulated BAI for half of the period
    # This is crucial in terms of correct harvesting and mortality estimates

    # Simulate BAI for halfPeriod
    initial_df <- BAI_prediction_halfPeriod(df_fit = data_BAI,
                                            df_predict = initial_df,
                                            site_vars = site_vars,
                                            rf_mtry = rf_mtry,
                                            species_n_threshold = species_n_threshold,
                                            include_climate = include_climate,
                                            measurement_thresholds = measurement_thresholds,
                                            area_correction = area_correction)

    # Next, we simulate height and crownHeight based on half period attributes

    # Calculate tree heights - half Period
    initial_df <- height_prediction_halfPeriod(df_fit = data_height, df_predict = initial_df,
                                               species_n_threshold = species_n_threshold,
                                               height_model = height_model,
                                               BRNN_neurons = BRNN_neurons,
                                               height_pred_level = height_pred_level)

    # Calculate tree crownHeights half Period
    if (sim_crownHeight == TRUE){

      initial_df <- crownHeight_prediction_halfPeriod(df_fit = data_crownHeight,
                                                      df_predict = initial_df,
                                                      site_vars = site_vars,
                                                      crownHeight_model = crownHeight_model,
                                                      BRNN_neurons = BRNN_neurons,
                                                      species_n_threshold = species_n_threshold)
    } else {

      initial_df$crownHeight_mid <- NA

    }

    # Calculate Volume
    if (volume_calculation == "form_factors"){

      initial_df <- vol_form_factors_halfPeriod(df = initial_df, form_factors = form_factors,
                                                form_factors_level = form_factors_level,
                                                uniform_form_factor = uniform_form_factor)

    } else if (volume_calculation == "volume_functions"){

      initial_df$volume_mid <- NA

      if (is.null(data_volF_param)){

        stop("data_volF_param is not provided")

      }

      initial_df <- V_general_halfPeriod(df = initial_df, data_volF_param)

    } else if (volume_calculation == "tariffs"){

      if (is.null(data_tariffs)){

        stop("data_tariffs is not supplied!")
      }

      initial_df$volume_mid <- NA

      initial_df <- vol_tariffs_halfPeriod(df = initial_df, data_tariffs = data_tariffs)

    } else if(volume_calculation == "slo_2p_volume_functions"){

      initial_df$volume_mid <- NA

      if (merchantable_whole_tree == "merchantable"){
        initial_df <- volume_merchantable_halfPeriod(df = initial_df)
      } else if (merchantable_whole_tree == "whole_tree"){
        initial_df <- volume_whole_tree_halfPeriod(df = initial_df)
      }

    } else {

      stop("Please define volume calculations: form_factors, volume_functions or tariffs")

    }


    # Ingrowth is now added, BAL and stand variables have changed. We therefore
    # update all variables
    initial_df <- calculate_BAL_halfPeriod(df = initial_df)
    initial_df <- add_stand_variables_halfPeriod(df = initial_df)

    # Save results

    if (sim_export_mode == TRUE){

      list_results[[sim]] <- initial_df

    }

    if (export_csv == TRUE){

      write.csv(initial_df, paste0("output_step_", sim-1, ".csv"), row.names = FALSE)

    }

    setTxtProgressBar(pb,sim) # progress bar

  }

  # remove data climate to save space
  rm(data_climate)

  # Select columns for the output
  final_calculations <- dplyr::select(do.call(bind_rows, list_results),
         "plotID", "treeID", "species", "speciesGroup", "year", "code",
         "weight", "p_weight", "weight_mid", "p_weight_mid",
         "height", "p_height", "height_mid","p_height_mid",
         "crownHeight", "p_crownHeight", "crownHeight_mid", "p_crownHeight_mid",
         "BA", "p_BA", "BA_mid", "p_BA_mid",
         "volume", "p_volume", "volume_mid", "p_volume_mid",
         "BAI","BAI_mid",
         "stand_BA", "stand_n", "BAL",
         site_vars)

  final_ouputs <- list(

    sim_results = final_calculations,

    height_eval = h_predictions$data_height_eval,
    crownHeight_eval = Crown_h_predictions$eval_crownHeight,
    mortality_eval = eval_mortality_output,
    ingrowth_eval = eval_ingrowth_output,
    BAI_eval = eval_BAI_output,

    height_model_species = h_predictions$model_species,
    height_model_speciesGroups = h_predictions$model_speciesGroups,
    crownHeight_model_species = Crown_h_predictions$model_species,
    crownHeight_model_speciesGroups = Crown_h_predictions$model_speciesGroups,
    mortality_model = mortality_output_model,
    BAI_model_species = BAI_outputs_model_species,
    BAI_model_speciesGroups = BAI_outputs_model_groups,
    max_size = max_size

  )

  list_n_elements <- length(final_ouputs)
  n_el <- 1

  # append the ingrowth models
  if (sim_ingrowth == TRUE){
    for (i in 1:length(ing_codes)){

	    final_ouputs[[list_n_elements + n_el]] <- get(paste0("ing_model_output_", ing_codes[i]))
      names(final_ouputs)[list_n_elements + n_el] <- paste0("ingrowth_model_",ing_codes[i])
	    n_el <- n_el + 1

    }
  }

  close(pb)

  return(final_ouputs)

}

