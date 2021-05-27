#' crownHeight_prediction
#'
#' Models to predict crown height
#'
#' @keywords internal
#'
crownHeight_prediction <- function(df_fit,  df_predict,
                                   site_vars = site_vars,
                                   species_n_threshold = 100,
                                   k = 10, eval_model_crownHeight = TRUE,
                                   crownHeight_model = "lm",
                                   BRNN_neurons = 3,
                                   blocked_cv = TRUE){

  # Define global variables
  plotID <- NULL
  treeID <- NULL
  year <- NULL
  BA <- NULL
  crownHeight <- NULL
  height <- NULL
  p_height <- NULL
  key_temp <- NULL
  key <- NULL
  speciesGroup <- NULL
  code <- NULL
  species <- NULL
  crownHeight_pred <- NULL

  df_fit <- mutate(df_fit, key = paste0(plotID, "_", treeID, "_", year))
  df_predict <- mutate(df_predict, key = paste0(plotID, "_", treeID, "_", year))

  # df <- rbind(df_fit, df_predict)
  # df <- df[!duplicated(df$key), ]

  ###################
  # crownHeight prediction
  ###################

  # Remove trees without BA
  no_BA <- filter(df_predict, is.na(BA))

  ##################### Loop for all traiff groups ########################
  # We define here strategy: if there are enough measurements, we calculate crownHeights for species
  # If not, we use grouping ID
  crownHeight_data <-  dplyr::filter(df_fit, !is.na(crownHeight))
  crownHeight_data <-  table(crownHeight_data$species)

  # Here we define instances that will be prediced using group variable
  data_below_threshold <- droplevels(df_fit[df_fit$species %in% names(crownHeight_data)[crownHeight_data < species_n_threshold],,drop=FALSE])
  uniq_tSk <- unique(data_below_threshold$speciesGroup)

  # Here we define species that will be predicted as individual species
  data_above_threshold <- droplevels(df_fit[df_fit$species %in% names(crownHeight_data)[crownHeight_data >= species_n_threshold],,drop=FALSE])
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

    dv_temporal_predict <- filter(dv_temporal_predict, !is.na(BA)) # Ta step sem dal ven, ker ni potreben

    temp_df <- select(dv_temporal_fit, crownHeight, height, all_of(site_vars))

    ########
    # eval #
    ########

    if (eval_model_crownHeight == TRUE){

      if (blocked_cv == FALSE){

        sampled_sequence <- sample(nrow(temp_df))

        df_fit_cv <- dv_temporal_fit[sampled_sequence, ]
        df_predict_cv <- temp_df[sampled_sequence, ]

      } else {

        df_fit_cv <- dv_temporal_fit
        df_predict_cv <- temp_df

      }

    foldi <- seq(1:k)
    folds <- cut(seq(1, nrow(temp_df)), breaks = k, labels = FALSE)

    for (m in 1:k){

      testIndexes <- which(folds == m, arr.ind = TRUE)

      test <- df_fit_cv[testIndexes, ]
      train <- df_predict_cv[-testIndexes, ]

      if (crownHeight_model == "lm"){

        ch_model <- lm(crownHeight ~ ., data = train)

      } else if(crownHeight_model == "brnn") {

        capture.output(ch_model <- brnn(crownHeight ~ ., data = train, neurons = BRNN_neurons))

      }

      test$crownHeight_pred <- predict(ch_model, test)

      eval_list[[m_holder]] <- test
      m_holder = m_holder + 1
      }

    }

    ###################
    # prediction part #
    ###################

    if (crownHeight_model == "lm"){

      mod1 <- lm(crownHeight ~ ., data = temp_df)

    } else if(crownHeight_model == "brnn") {

      capture.output(mod1 <- brnn(crownHeight ~ ., data = temp_df, neurons = BRNN_neurons))

    }

    dv_temporal_predict$crownHeight <- predict(mod1, dv_temporal_predict)

    # in predict.brnn na.pass is not working, so I do workaround here
    dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))
    dv_temporal_predict_noNA <- filter(dv_temporal_predict, !is.na(p_height))
    dv_temporal_predict_yesNA <- filter(dv_temporal_predict, is.na(p_height))
    dv_temporal_predict_noNA$p_crownHeight <- predict(mod1, newdata = rename(dv_temporal_predict_noNA,
                                                                             "height_X" = "height",
                                                                             "height" = "p_height"))

    dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dv_temporal_predict_yesNA)
    dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
    dv_temporal_predict$key_temp <- NULL

    list_predictions[[p]] <- dv_temporal_predict
    p = p + 1

  }

  DF_predictions_species <- do.call(rbind, list_predictions)

  DF_eval_species <- do.call(rbind, eval_list)

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
    dv_temporal_predict <- filter(dv_temporal_predict, !is.na(BA))

    temp_df <- select(dv_temporal_fit, crownHeight, height, all_of(site_vars))

    ########
    # eval #
    ########
    if (eval_model_crownHeight == TRUE){

      if (blocked_cv == FALSE){

        sampled_sequence <- sample(nrow(temp_df))

        df_fit_cv <- dv_temporal_fit[sampled_sequence, ]
        df_predict_cv <- temp_df[sampled_sequence, ]

      } else {

        df_fit_cv <- dv_temporal_fit
        df_predict_cv <- temp_df

      }

    foldi <- seq(1:k)
    folds <- cut(seq(1, nrow(temp_df)), breaks = k, labels = FALSE)

    for (m in 1:k){

      testIndexes <- which(folds == m, arr.ind = TRUE)

      test <- df_fit_cv[testIndexes, ]
      train <- df_predict_cv[-testIndexes, ]

      if (crownHeight_model == "lm"){

        ch_model <- lm(crownHeight ~ ., data = train)

      } else if(crownHeight_model == "brnn") {

        capture.output(ch_model <- brnn(crownHeight ~ ., data = train, neurons = BRNN_neurons))

      }

      test$crownHeight_pred <- predict(ch_model, test)

      eval_list[[m_holder]] <- test
      m_holder = m_holder + 1

    }
    }
    ###################
    # prediction part #
    ###################


    if (crownHeight_model == "lm"){

      mod1 <- lm(crownHeight ~ ., data = temp_df)

    } else if(crownHeight_model == "brnn") {

      capture.output(mod1 <- brnn(crownHeight ~ ., data = temp_df, neurons = BRNN_neurons))

    }

    dv_temporal_predict$crownHeight <- predict(mod1, dv_temporal_predict)

    # in predict.brnn na.pass is not working, so I do workaround here
    dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))
    dv_temporal_predict_noNA <- filter(dv_temporal_predict, !is.na(p_height))
    dv_temporal_predict_yesNA <- filter(dv_temporal_predict, is.na(p_height))
    dv_temporal_predict_noNA$p_crownHeight <- predict(mod1, newdata = rename(dv_temporal_predict_noNA,
                                                                             "height_X" = "height",
                                                                             "height" = "p_height"))

    dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dv_temporal_predict_yesNA)
    dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
    dv_temporal_predict$key_temp <- NULL

    list_predictions[[p]] <- dv_temporal_predict
    p = p + 1

  }

  DF_predictions_sGroups <- do.call(rbind, list_predictions)
  DF_predictions_sGroups <-filter(DF_predictions_sGroups, !(key %in% DF_predictions_species$key))

  DF_eval_sGroups <- do.call(rbind, eval_list)

  if (eval_model_crownHeight == TRUE){
    DF_eval_sGroups <-filter(DF_eval_sGroups, !(key %in% DF_eval_species$key))
    data_eval_crownHeight_predictions <- rbind(
      DF_eval_sGroups,
      DF_eval_species)

    data_eval_crownHeight_predictions <- select(data_eval_crownHeight_predictions, plotID, treeID, year, speciesGroup, code,
                                                species, crownHeight, crownHeight_pred)

  } else {

    data_eval_crownHeight_predictions <- "the argument set_eval_crownHeight is set to FALSE"

  }

  # Now add crownHeight data to our data (merge)
  data_crownHeight_predictions <- rbind(
    DF_predictions_sGroups,
    DF_predictions_species)

  data_crownHeight_predictions <- rbind(data_crownHeight_predictions, no_BA)

  # correction: Negative crownHeight predictions are set to 1 meter (could also be the minimum measured)
  data_crownHeight_predictions <- mutate(data_crownHeight_predictions, crownHeight = ifelse(crownHeight < 0, 1, crownHeight))

  data_crownHeight_predictions$key <- NULL

  ##########################################

  final_output_list <- list(

    predicted_crownHeight = data_crownHeight_predictions,
    eval_crownHeight = data_eval_crownHeight_predictions

  )

  return(final_output_list)

}

