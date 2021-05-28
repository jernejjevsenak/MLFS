#' height_prediction
#'
#' Height model
#'
#' @param df_fit data frame with data being used for model fit
#' @keywords internal
#'

height_prediction <- function(df_fit,  df_predict,
                              species_n_threshold = 100,
                              height_model = "naslund",
                              BRNN_neurons = 3,
                              height_pred_level = 0,
                              eval_model_height = TRUE,
                              blocked_cv = TRUE, k = 10
                              ){

  # Define global variables
  plotID <- NULL
  treeID <- NULL
  year <- NULL
  BA <- NULL
  height <- NULL
  p_BA <- NULL
  key_temp <- NULL
  key <- NULL
  speciesGroup <- NULL
  code <- NULL
  species <- NULL
  h_pred <- NULL

  df_fit <- mutate(df_fit, key = paste0(plotID, "_", treeID, "_", year))
  df_predict <- mutate(df_predict, key = paste0(plotID, "_", treeID, "_", year))

  ###################
  # Height prediction
  ###################

  # Remove trees withour BA
  # no_BA <- filter(df_predict, is.na(BA))

  ##################### Loop for all traiff groups ########################
  # We define here strategy: if there are enough measurements, we calculate heights for species
  # If not, we use grouping ID
  height_data <-  dplyr::filter(df_fit, !is.na(height))
  height_data <-  table(height_data$species)

  # Here we define instances that will be prediced using group variable
  data_below_threshold <- droplevels(df_fit[df_fit$species %in% names(height_data)[height_data < species_n_threshold],,drop=FALSE])
  uniq_tSk <- unique(data_below_threshold$speciesGroup)

  # Here we define species that will be predicted as individual species
  data_above_threshold <- droplevels(df_fit[df_fit$species %in% names(height_data)[height_data >= species_n_threshold],,drop=FALSE])
  unique_dv <- unique(data_above_threshold$species)

  ######################################################################
  ######################################################################

  list_predictions <- list()
  p = 1

  eval_list <- list()
  m_holder = 1

  for (M in unique_dv){

    ###################### izbira dv ########################################
    dv_temporal_fit <- subset(df_fit, subset = df_fit$species %in% M)
    # BA can not be missing!
    dv_temporal_fit <- filter(dv_temporal_fit, !is.na(BA))


    dv_temporal_predict <- subset(df_predict, subset = df_predict$species %in% M)
    # dv_temporal_predict <- filter(dv_temporal_predict, !is.na(BA))

    if (height_model == "brnn"){

      if (height_pred_level == 1){

        capture.output(mod1 <- brnn(height ~ BA + plotID, data = dv_temporal_fit, neurons = BRNN_neurons))

      } else {

        capture.output(mod1 <- brnn(height ~ BA, data = dv_temporal_fit, neurons = BRNN_neurons))

      }

    } else {

      mod1 <- fithd(dv_temporal_fit$BA,
                    dv_temporal_fit$height,
                    plot = dv_temporal_fit$plotID,
                    modelName = height_model, control=list(maxIter = 1000, msmaxIter = 1000))

    }


    ##############
    # eval_step #
    #############

    if (eval_model_height == TRUE){


    if (height_model == "brnn"){

      if (blocked_cv == FALSE){

        sampled_sequence <- sample(nrow(dv_temporal_fit))

        df_fit_cv <- dv_temporal_fit[sampled_sequence, ]

      } else {

        df_fit_cv <- dv_temporal_fit

      }

      foldi <- seq(1:k)
      folds <- cut(seq(1, nrow(df_fit_cv)), breaks = k, labels = FALSE)

      for (m in 1:k){

        testIndexes <- which(folds == m, arr.ind = TRUE)

        test <- df_fit_cv[testIndexes, ]
        train <- df_fit_cv[-testIndexes, ]

        if(height_pred_level == 1) {

          capture.output(mod_brnn <- brnn(height ~ BA + plotID, data = train, neurons = BRNN_neurons))

        } else if (height_pred_level == 0) {

          capture.output(mod_brnn <- brnn(height ~ BA, data = train, neurons = BRNN_neurons))

        }

        test$h_pred <- predict(mod_brnn, newdata = test)

        eval_list[[m_holder]] <- test
        m_holder = m_holder + 1

      }

    } else {

      eval_step <- dv_temporal_fit

      eval_step$h_pred <- predict(mod1, newdata = rename(dv_temporal_predict, "d" = "BA", "plot" = "plotID"), level = height_pred_level)

      eval_list[[p]] <- eval_step

    }

  }

    ###################
    # prediction step #
    ###################

    if (height_model == "brnn"){

      # in predict.brnn na.pass is not working, so I do workaround here

      # 1 predictions for BA
      dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))
      dv_temporal_predict_noNA <- filter(dv_temporal_predict, !is.na(BA))
      dv_temporal_predict_yesNA <- filter(dv_temporal_predict, is.na(BA))
      dv_temporal_predict_noNA$height <- predict(mod1, newdata =dv_temporal_predict_noNA)
      dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dv_temporal_predict_yesNA)
      dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
      dv_temporal_predict$key_temp <- NULL

      # 2 predictions for p_BA
      dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))
      dv_temporal_predict_noNA <- filter(dv_temporal_predict, !is.na(p_BA))
      dv_temporal_predict_yesNA <- filter(dv_temporal_predict, is.na(p_BA))
      dv_temporal_predict_noNA$p_height <- predict(mod1,
                                                   newdata = rename(dv_temporal_predict_noNA,
                                                                          "BA_temp" = "BA",
                                                                          "BA" = "p_BA"))

      dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dv_temporal_predict_yesNA)
      dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
      dv_temporal_predict$key_temp <- NULL

    } else {

      dv_temporal_predict$height <- predict(mod1, newdata = rename(dv_temporal_predict, "d" = "BA", "plot" = "plotID"), level = height_pred_level)
      dv_temporal_predict$p_height <- predict(mod1, newdata = rename(dv_temporal_predict,
                                                                       "d" = "p_BA", "plot" = "plotID"), level = height_pred_level, na.action = na.pass)
    }

    list_predictions[[p]] <- dv_temporal_predict
    p = p + 1

  }

  DF_predictions_species <- do.call(rbind, list_predictions)
  DF_eval_species <-  do.call(rbind, eval_list)

####################################################
# Now we repeat the process for minor tree species #
####################################################

  list_predictions <- list()
  p = 1

  eval_list <- list()
  m_holder = 1

  for (M in uniq_tSk){

    ###################### select dv #############################
    dv_temporal_fit <- subset(df_fit, subset = df_fit$speciesGroup %in% M)

    # BA can not be missing!
    dv_temporal_fit <- filter(dv_temporal_fit, !is.na(BA))

    dv_temporal_predict <- subset(df_predict, subset = df_predict$speciesGroup %in% M)
    # dv_temporal_predict <- filter(dv_temporal_predict, !is.na(BA))


    if (height_model == "brnn"){

      if (height_pred_level == 1){

        capture.output(mod1 <- brnn(height ~ BA + plotID, data = dv_temporal_fit, neurons = BRNN_neurons))

      } else {

        capture.output(mod1 <- brnn(height ~ BA, data = dv_temporal_fit, neurons = BRNN_neurons))

      }

    } else {

      mod1 <- fithd(dv_temporal_fit$BA,
                    dv_temporal_fit$height,
                    plot = dv_temporal_fit$plotID,
                    modelName = height_model, control=list(maxIter = 1000, msmaxIter = 1000))

    }



    #############
    # eval_step #
    #############

    if (eval_model_height == TRUE){


      if (height_model == "brnn"){

        if (blocked_cv == FALSE){

          sampled_sequence <- sample(nrow(dv_temporal_fit))

          df_fit_cv <- dv_temporal_fit[sampled_sequence, ]

        } else {

          df_fit_cv <- dv_temporal_fit

        }

        foldi <- seq(1:k)
        folds <- cut(seq(1, nrow(df_fit_cv)), breaks = k, labels = FALSE)

        for (m in 1:k){

          testIndexes <- which(folds == m, arr.ind = TRUE)

          test <- df_fit_cv[testIndexes, ]
          train <- df_fit_cv[-testIndexes, ]

          if(height_pred_level == 1) {

            capture.output(mod_brnn <- brnn(height ~ BA + plotID, data = train, neurons = BRNN_neurons))

          } else if (height_pred_level == 0) {

            capture.output(mod_brnn <- brnn(height ~ BA, data = train, neurons = BRNN_neurons))

          }

          test$h_pred <- predict(mod_brnn, newdata = test)

          eval_list[[m_holder]] <- test
          m_holder = m_holder + 1

        }

      } else {

        eval_step <- dv_temporal_fit

        eval_step$h_pred <- predict(mod1, newdata = rename(dv_temporal_predict, "d" = "BA", "plot" = "plotID"), level = height_pred_level)

        eval_list[[p]] <- eval_step

      }

    }

    ###################
    # prediction step #
    ###################

    if (height_model == "brnn"){

      # in predict.brnn na.pass is not working, so I do workaround here

      # 1 predictions for BA
      dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))
      dv_temporal_predict_noNA <- filter(dv_temporal_predict, !is.na(BA))
      dv_temporal_predict_yesNA <- filter(dv_temporal_predict, is.na(BA))
      dv_temporal_predict_noNA$height <- predict(mod1, newdata =dv_temporal_predict_noNA)
      dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dv_temporal_predict_yesNA)
      dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
      dv_temporal_predict$key_temp <- NULL

      # 2 predictions for p_BA
      dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))
      dv_temporal_predict_noNA <- filter(dv_temporal_predict, !is.na(p_BA))
      dv_temporal_predict_yesNA <- filter(dv_temporal_predict, is.na(p_BA))
      dv_temporal_predict_noNA$p_height <- predict(mod1,
                                                   newdata = rename(dv_temporal_predict_noNA,
                                                                    "BA_temp" = "BA",
                                                                    "BA" = "p_BA"))

      dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dv_temporal_predict_yesNA)
      dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
      dv_temporal_predict$key_temp <- NULL

    } else {

      dv_temporal_predict$height <- predict(mod1, newdata = rename(dv_temporal_predict, "d" = "BA", "plot" = "plotID"), level = height_pred_level)
      dv_temporal_predict$p_height <- predict(mod1, newdata = rename(dv_temporal_predict,
                                                                       "d" = "p_BA", "plot" = "plotID"), level = height_pred_level, na.action = na.pass)
      }

    list_predictions[[p]] <- dv_temporal_predict
    p = p + 1

  }

  DF_predictions_sGroups <- do.call(rbind, list_predictions)
  DF_predictions_sGroups <-filter(DF_predictions_sGroups, !(key %in% DF_predictions_species$key))

  if (eval_model_height == TRUE){

  DF_eval_sGroups <- do.call(rbind, eval_list)
  DF_eval_sGroups <-filter(DF_eval_sGroups, !(key %in% DF_eval_species$key))

  data_height_eval <- rbind(
    DF_eval_sGroups,
    DF_eval_species)

  data_height_eval <- select(data_height_eval, plotID, treeID, year, speciesGroup, code, species, height, BA, h_pred)

  } else {

    data_height_eval <- "the argument set_eval_height is set to FALSE"

  }

  # Now add height data to our data (merge)
  data_height_predictions <- rbind(
    DF_predictions_sGroups,
    DF_predictions_species)

  # data_height_predictions <- rbind(data_height_predictions, no_BA)

  # correction: Negative height predictions are set to 1 meter (could also be the minimum measured)
  data_height_predictions <- mutate(data_height_predictions, height = ifelse(height < 0, 1, height))

  data_height_predictions$key <- NULL

  #####################################################

  final_output_list <- list(

    data_height_predictions = data_height_predictions,
    data_height_eval = data_height_eval

  )

  return(final_output_list)

}

