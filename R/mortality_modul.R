#' predict_mortality
#'
#' This sub model first fits a binary model to derive the effects of individual
#' tree, site and climate variables on mortality; and afterwards predict the
#' probability of dying for each tree from df_predict
#'
#' @param df_fit a data frame with individual tree data and site descriptors
#' where code is used to specify a status of each tree
#' @param df_predict data frame which will be used for mortality predictions
#' @param df_climate data frame with monthly climate data
#' @param mortality_share a value defining the proportion of the volume which is
#' to be the subject of mortality
#' @param mortality_share_type character, it can be 'volume' or 'n_trees'. If
#' 'volume' then the mortality share relates to total standing volume, if
#' 'n_trees' then mortality share relates to the total number of standing trees
#' @param include_climate logical, should climate variables be included as
#' predictors
#' @param site_vars a character vector of variable names which are used as site
#' descriptors
#' @param select_months_climate vector of subset months to be considered.
#' Default is c(1,12), which uses all months.
#' @param mortality_model logical, should the mortality model be evaluated
#' and returned as the output
#' @param nb_laplace value used for Laplace smoothing (additive smoothing) in
#' naive Bayes algorithm. Defaults to 0 (no Laplace smoothing).
#' @param sim_crownHeight logical, should crown heights be considered as a
#' predictor variable? If TRUE, a crownHeight column is expected in data_NFI
#' @param k the number of folds to be used in the k fold cross-validation
#' @param eval_model_mortality logical, should the mortality model be evaluated
#' and returned as the output
#' @param blocked_cv logical, should the blocked cross-validation be used in the
#' evaluation phase?
#' @param sim_mortality logical, should mortality be simulated?
#' @param sim_step_years the simulation step in years
#' @param rf_mtry number of variables randomly sampled as candidates at each
#' split of a random forest model. If NULL, default settings are applied.
#' @param df_max_size a data frame with the maximum BA values for each species.
#' If a tree exceeds this value, it dies.
#' @param ingrowth_codes numeric value or a vector of codes which refer to
#' ingrowth trees
#' @param include_mortality_BAI logical, should basal area increments (BAI) be
#' used as independent variable for predicting individual tree morality?
#' @param intermediate_print logical, if TRUE intermediate steps will be printed
#' while the mortality sub model is running
#' @param use_max_size_threshold logical - should the principle of maxium size
#' be applied?
#' @param mortality_bias_adjusted Logical (length-one). If `TRUE` (default),
#'   applies a simple bias fix so large trees aren’t over-removed. The frequency
#'   of adjustment is controlled by the `bias_adj_factor` argument. If `FALSE`,
#'   predicted probabilities are left unchanged.
#' @param bias_adj_factor Integer (>= 2). Controls how sparsely you reduce
#' death probabilities among the top-ranked trees. Starting from the 3rd row,
#' every `bias_adj_factor`-th tree has its probability set to zero—so `2` keeps
#' every second high-risk tree alive, `3` every third, and so on.
#'
#' @return a list with three elements:
#' \enumerate{
#'  \item $predicted_mortality - a data frame with updated tree status (code) based on the predicted mortality
#'  \item $eval_mortality - a data frame with predicted and observed probabilities of dying for all individual trees, or character string indicating that mortality sub-model was not evaluated
#'  \item $model_output - the output model for mortality
#' }
#'
#' @examples
#' data("data_v4")
#' data("data_mortality")
#' data("max_size_data")
#'
#' mortality_outputs <- predict_mortality(
#'  df_fit = data_mortality,
#'  df_predict = data_v4,
#'  mortality_share_type = 'volume',
#'  df_climate = data_climate,
#'  site_vars = c("slope", "elevation", "northness", "siteIndex"),
#'  sim_mortality = TRUE,
#'  mortality_model = 'naiveBayes',
#'  nb_laplace = 0,
#'  sim_crownHeight = TRUE,
#'  mortality_share = 0.02,
#'  include_climate = TRUE,
#'  select_months_climate = c(6,7,8),
#'  eval_model_mortality = TRUE,
#'  k = 10, blocked_cv = TRUE,
#'  sim_step_years = 6,
#'  df_max_size = max_size_data,
#'  ingrowth_codes = c(3,15),
#'  include_mortality_BAI = TRUE)
#'
#'  df_predicted <- mortality_outputs$predicted_mortality
#'  df_evaluation <- mortality_outputs$eval_mortality
#'
#'  # confusion matrix
#'  table(df_evaluation$mortality, round(df_evaluation$mortality_pred, 0))
#'

predict_mortality <- function(df_fit, df_predict, df_climate, mortality_share = NA,
                              mortality_share_type = "volume",
                              include_climate, site_vars, select_months_climate = c(6,8),
                              mortality_model = "rf", nb_laplace = 0, sim_crownHeight = FALSE,
                              k = 10, eval_model_mortality = TRUE, blocked_cv = TRUE,
                              sim_mortality = TRUE, sim_step_years = 5, rf_mtry = NULL,
                              df_max_size = NULL, ingrowth_codes = 3, include_mortality_BAI = TRUE,
                              intermediate_print = FALSE, use_max_size_threshold = FALSE,
                              mortality_bias_adjusted = TRUE, bias_adj_factor = 2


                              ){

# Define global variables
year <- NULL
month <- NULL
plotID <- NULL
p_sum <- NULL
t_avg <- NULL
code <- NULL
treeID <- NULL
speciesGroup <- NULL
species <- NULL
mortality <- NULL
mortality_pred <- NULL
p_mortality <- NULL
BA <- NULL
BAL <- NULL
height <- NULL
crownHeight <- NULL
stand_BA <- NULL
stand_n <- NULL
col_sum <- NULL
vol_ha_mid <- NULL
BA_mid <- NULL
BA_max <- NULL
ranger <- NULL
BAI <- NULL
BAI_mid <- NULL
n_obs <- NULL
BA_bin <- NULL
n_tgt <- NULL


if (sim_mortality == TRUE){

  if (include_climate == TRUE){

    initial_colnames <- colnames(df_predict)

    max_year_predict <- max(df_predict$year)
    min_year_predict <- max_year_predict + sim_step_years

    climate_predict <- dplyr::filter(df_climate, year %in% seq(min_year_predict, max_year_predict)) %>%
      dplyr::filter(month %in% select_months_climate) %>%
      group_by(plotID) %>% summarise(p_sum = sum(p_sum), t_avg = mean(t_avg))

    df_predict$t_avg <- NULL
    df_predict$p_sum <- NULL

    df_predict <- merge(df_predict, climate_predict, by = "plotID")

    df_predict <- dplyr::select(df_predict, all_of(initial_colnames))

    site_vars <- c(site_vars, "p_sum", "t_avg")

  }

  df_predict <- dplyr::filter(df_predict, code %in% c(0,ingrowth_codes))

  if (sim_crownHeight == TRUE){

    df_fit <- mutate(df_fit, mortality = ifelse(code == 2, 1, 0),
                     BA_mid = BA,
                     BAL_mid = BAL,
                     BAI_mid = BAI,
                     height_mid = height,
                     crownHeight_mid = crownHeight,
                     stand_BA_mid = stand_BA,
                     stand_n_mid = stand_n)

  } else {

    df_fit <- mutate(df_fit, mortality = ifelse(code == 2, 1, 0),
                     BA_mid = BA,
                     BAL_mid = BAL,
                     height_mid = height,
                     # crownHeight_mid = crownHeight,
                     stand_BA_mid = stand_BA,
                     stand_n_mid = stand_n)

  }

  # YOU CAN MANUALLY SET THE MORTALITY SHARE
  if (is.na(mortality_share)){

    mortality_share <- sum(df_fit$code == 2)/nrow(df_fit)

    if (intermediate_print == TRUE){

      message(paste0("Estimated moratlity share is ", round(mortality_share, 2)))

    }

  }

  if (sim_crownHeight == TRUE){

    if (include_mortality_BAI == TRUE){

      formula <- as.formula(paste0("mortality ~ BA_mid + BAI_mid + height_mid + crownHeight_mid + BAL_mid + stand_BA_mid + stand_n_mid + speciesGroup +",
                                   paste(all_of(site_vars), collapse = "+")))

      df_fit <- filter(df_fit, !is.na(BAI_mid))

    } else {

      formula <- as.formula(paste0("mortality ~ BA_mid + height_mid + crownHeight_mid + BAL_mid + stand_BA_mid + stand_n_mid + speciesGroup +",
                                   paste(all_of(site_vars), collapse = "+")))

    }

  } else {

    if (include_mortality_BAI == TRUE){

    formula <- as.formula(paste0("mortality ~ BA_mid + BAI_mid + height_mid + BAL_mid + stand_BA_mid + stand_n_mid + speciesGroup +",
                                 paste(all_of(site_vars), collapse = "+")))

    df_fit <- filter(df_fit, !is.na(BAI_mid))

    summary(df_fit)

    } else {

      formula <- as.formula(paste0("mortality ~ BA_mid + height_mid + BAL_mid + stand_BA_mid + stand_n_mid + speciesGroup +",
                                   paste(all_of(site_vars), collapse = "+")))


    }

  }


  ##############
  # Eval phase #
  ##############

  if (eval_model_mortality == TRUE){

    foldi <- seq(1:k)
    folds <- cut(seq(1, nrow(df_fit)), breaks = k, labels = FALSE)

    eval_list <- list()

    if (blocked_cv == FALSE){

      df_fit_cv <- df_fit[sample(nrow(df_fit)), ]

      } else {

      df_fit_cv <- df_fit
    }

    for (m in 1:k){

      testIndexes <- which(folds == m, arr.ind = TRUE)

      test <- df_fit_cv[testIndexes, ]
      train <- df_fit_cv[-testIndexes, ]

      if (mortality_model == "glm"){

        model_mortality <- glm(formula, data = train, family = "binomial")
        test$mortality_pred <- predict(model_mortality, test, type="response")


      } else if (mortality_model == "rf"){

        if (is.null(rf_mtry)){

          model_mortality <- ranger(formula, data = train)

        } else {

          model_mortality <- ranger(formula, data = train, mtry = rf_mtry)

        }

        test$mortality_pred <- predict(model_mortality, test, type="response")$predictions

      } else if (mortality_model == "naiveBayes") {

        train$mortality <- factor(train$mortality)
        model_mortality <- naive_bayes(formula, data = train, laplace = nb_laplace)
        suppressWarnings(test$mortality_pred <- predict(model_mortality, test, type="prob")[,2])

      } else {

        stop("mortality_model should be 'glm', 'naiveBayes' or 'rf'")

      }

      eval_list[[m]] <- test

    }

    df_eval_mortality <- do.call(rbind, eval_list)

    df_eval_mortality <- dplyr::select(df_eval_mortality, plotID, treeID, year, speciesGroup, code,
                                species, mortality,  mortality_pred)
  } else {

    df_eval_mortality <- "the argument set_eval_mortality is set to FALSE"

  }






  ####################
  # Prediction phase #
  ####################

  if (mortality_model == "glm"){

    model_mortality <- glm(formula, data = df_fit, family = "binomial")
    df_predict$p_mortality <- predict(model_mortality, df_predict, type="response")

  } else if (mortality_model == "rf"){

    if (is.null(rf_mtry)){

      model_mortality <- ranger(formula, data = df_fit)

    } else {

      model_mortality <- ranger(formula, data = df_fit, mtry = rf_mtry)

    }

    df_predict$p_mortality <- predict(model_mortality, df_predict, type="response")$predictions

  } else  if (mortality_model == "naiveBayes"){

    df_fit$mortality <- factor(df_fit$mortality)
    model_mortality <- naive_bayes(formula, data = df_fit, laplace = nb_laplace)
    suppressWarnings(df_predict$p_mortality <- predict(model_mortality, df_predict, type="prob")[,2])

  } else {

    stop("mortality_model should be 'glm', 'naiveBayes' or 'rf'")

  }




















  if (use_max_size_threshold == TRUE){

    df_predict <- merge(df_predict, df_max_size, by = 'species', all.x = TRUE)

    df_predict <- dplyr::mutate(df_predict, p_mortality = ifelse(BA_mid > BA_max, 1, p_mortality),
                                BA_max = NULL
    ) %>%
      arrange(-p_mortality)

  } else {

    df_predict <- df_predict %>%
      arrange(-p_mortality)

  }





  #################
  # Balance data #
  ################

  # Mortality model could bias mortality predictions towards trees with larger DBH
  if (mortality_bias_adjusted == TRUE){

    # assume numeric p_mortality and descending order already applied
    p <- df_predict$p_mortality

    # pick rows 3, 6, 9, ...
    idx <- seq(3, length(p), by = bias_adj_factor)

    # df_predict$p_mortality_adj <- p
    df_predict$p_mortality[idx] <- 0

    df_predict <- df_predict %>%
      arrange(-p_mortality)

  }



  if (mortality_share_type == "volume"){

    df_predict$vol_ha_mid <- df_predict$volume_mid * df_predict$weight_mid
    volume_total <- sum(df_predict$vol_ha_mid, na.rm = TRUE)
    df_predict <- mutate(df_predict,
                         col_sum = cumsum(replace_na(vol_ha_mid, 0)),
                         code = ifelse(col_sum < (volume_total * mortality_share), 2, 0))

    df_predict[, "year"] <- df_predict[, "year"] + sim_step_years

    df_predict$p_mortality <- NULL
    df_predict$col_sum <- NULL
    df_predict$vol_ha_mid <- NULL

  } else if (mortality_share_type == "n_trees"){

    cut_th <- round(nrow(df_predict) * mortality_share)
    df_predict[c(1:cut_th), "code"] <- 2
    df_predict[c((cut_th+1):nrow(df_predict)), "code"] <- 0

    df_predict[, "year"] <- df_predict[, "year"] + sim_step_years

  } else {

    stop(paste0("mortality_share_type should be 'n_trees' or 'volume' but instead it is ", mortality_share_type))

  }


  } else if (sim_mortality == FALSE){

    df_predict <- dplyr::filter(df_predict, code %in% c(0,3,15))
    df_predict[, "year"] <- df_predict[, "year"] + sim_step_years
    df_predict[ , "code"] <- 0
    df_eval_mortality <- paste0("sim_mortality is set to FALSE.",
    "Mortality is not simulated. eval_mortality is not available.")

    model_mortality <- paste0("sim_mortality is set to FALSE.",
                                "Mortality is not simulated. model_mortality is not available.")
    }

  df_predict <- arrange(df_predict, plotID, treeID)

  final_output_list <- list(

    predicted_mortality = df_predict,
    eval_mortality = df_eval_mortality,
    model_output = model_mortality

  )

  return(final_output_list)

}

