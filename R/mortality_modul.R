#' predict_mortality
#'
#' Mortality model
#'
#' @keywords internal
#'

predict_mortality <- function(df_fit, df_predict, df_climate, mortality_share = NA,
                              mortality_share_type = "volume",
                              include_climate, site_vars, select_months_climate = c(1:12),
                              mortality_model = "rf", nb_laplace = 0,
                              k = 10, eval_model_mortality = TRUE, blocked_cv = TRUE,
                              sim_mortality = TRUE, sim_step_years = 5, rf_mtry = NULL,
                              df_max_size = NULL){



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

  df_predict <- dplyr::filter(df_predict, code %in% c(0,3,15))
  df_fit <- mutate(df_fit, mortality = ifelse(code == 2, 1, 0),
                   BA_mid = BA,
                   BAL_mid = BAL,
                   height_mid = height,
                   crownHeight_mid = crownHeight,
                   stand_BA_mid = stand_BA,
                   stand_n_mid = stand_n)

  # YOU CAN MANUALLY SET THE MORTALITY SHARE
  if (is.na(mortality_share)){

    mortality_share <- sum(df_fit$code == 2)/nrow(df_fit)

    print(paste0("Estimated moratlity share is ", round(mortality_share, 2)))

  }

  formula <- as.formula(paste0("mortality ~ BA_mid + height_mid + crownHeight_mid + BAL_mid + stand_BA_mid + stand_n_mid + speciesGroup +",
                               paste(all_of(site_vars), collapse = "+")))

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

          model_mortality <- randomForest(formula, data = train)

        } else {

          model_mortality <- randomForest(formula, data = train, mtry = rf_mtry)

        }

        test$mortality_pred <- predict(model_mortality, test, type="response")

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

      model_mortality <- randomForest(formula, data = df_fit)

    } else {

      model_mortality <- randomForest(formula, data = df_fit, mtry = rf_mtry)

    }

    df_predict$p_mortality <- predict(model_mortality, df_predict, type="response")

  } else  if (mortality_model == "naiveBayes"){

    df_fit$mortality <- factor(df_fit$mortality)
    model_mortality <- naive_bayes(formula, data = df_fit, laplace = nb_laplace)
    suppressWarnings(df_predict$p_mortality <- predict(model_mortality, df_predict, type="prob")[,2])

  } else {

    stop("mortality_model should be 'glm', 'naiveBayes' or 'rf'")

  }

  df_predict <- merge(df_predict, df_max_size, by = 'species', all.x = TRUE)

  df_predict <- dplyr::mutate(df_predict, p_mortality = ifelse(BA_mid > BA_max, 1, p_mortality),
                              BA_max = NULL
                              ) %>%
    arrange(-p_mortality)

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

