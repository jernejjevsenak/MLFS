#' height_prediction_halfPeriod
#'
#' Height model
#'
#' @param df_fit data frame with data being used for model fit
#' @keywords internal
#'

height_prediction_halfPeriod <- function(df_fit,  df_predict,
                              species_n_threshold = 100,
                              height_model = "naslund",
                              BRNN_neurons = 3,
                              height_pred_level = 0
                              ){

  # Define global variables
  plotID <- NULL
  treeID <- NULL
  year <- NULL
  BA <- NULL
  BA_mid <- NULL
  height <- NULL
  height_mid <- NULL
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
  # no_BA <- dplyr::filter(df_predict, is.na(BA_mid))

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

  # check if there is any unique_speciesGroup missing
  if (any(!(unique(df_predict$speciesGroup)%in% uniq_tSk))){

    missing_spG <-  unique(df_predict$speciesGroup)[!(unique(df_predict$speciesGroup)%in% uniq_tSk)]
    unique_species_missing_spG <- unique(df_predict[df_predict$speciesGroup==missing_spG,"species"])

    # are these already in unique_dv?
    if (any(!(missing_spG %in% unique_dv))){

      # If not, append it
      uniq_tSk <- c(missing_spG, uniq_tSk)

    }
  }


  ######################################################################
  ######################################################################

  list_predictions <- list()
  p = 1

  for (M in unique_dv){

    ###################### izbira dv ########################################
    dv_temporal_fit <- subset(df_fit, subset = df_fit$species %in% M)


    # BA can not be missing!
    dv_temporal_fit <- dplyr::filter(dv_temporal_fit, !is.na(BA))

    # I change the name of BA into BA_mid, so that the predict function will select the right variable
    dv_temporal_fit <- rename(dv_temporal_fit, "BA_mid" = "BA")

    dv_temporal_predict <- dplyr::filter(df_predict, species %in% M)


    if (height_model == "brnn"){

      if (height_pred_level == 1){

        capture.output(mod1 <- brnn(height ~ BA_mid + plotID, data = dv_temporal_fit, neurons = BRNN_neurons))

      } else {

        capture.output(mod1 <- brnn(height ~ BA_mid, data = dv_temporal_fit, neurons = BRNN_neurons))

      }

    } else {

      mod1 <- lmfor::fithd(dv_temporal_fit$BA_mid,
                    dv_temporal_fit$height,
                    plot = dv_temporal_fit$plotID,
                    modelName = height_model, control=list(maxIter = 1000, msmaxIter = 1000))

    }

    ###################
    # prediction step #
    ###################

    if (height_model == "brnn"){

      # in predict.brnn na.pass is not working, so I do workaround here

      # 1 predictions for BA
      dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))

      dv_temporal_predict_noNA <- dplyr::filter(dv_temporal_predict, !is.na(BA_mid))
      dv_temporal_predict_yesNA <- dplyr::filter(dv_temporal_predict, is.na(BA_mid))
      dv_temporal_predict_noNA$height_mid <- predict(mod1, newdata =dv_temporal_predict_noNA)


      if (sum(is.na(dv_temporal_predict$BA)) > 0){

        dv_temporal_predict_yesNA$height_mid <- NA
        dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dplyr::select(dv_temporal_predict_yesNA, colnames(dv_temporal_predict_noNA)))

      } else {

        dv_temporal_predict <- dv_temporal_predict_noNA

      }

      dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
      dv_temporal_predict$key_temp <- NULL
      dv_temporal_predict$height_new <- NULL
      dv_temporal_predict$p_height_new <- NULL

    } else {

      dv_temporal_predict$height_mid <- predict(mod1, newdata = rename(dv_temporal_predict, "d" = "BA_mid", "plot" = "plotID"), level = height_pred_level, na.action = na.pass)

      dv_temporal_predict$height_new <- NULL
      dv_temporal_predict$p_height_new <- NULL

    }

    list_predictions[[p]] <- dv_temporal_predict
    p = p + 1

  }

  DF_predictions_species <- do.call(rbind, list_predictions)

####################################################
# Now we repeat the process for minor tree species #
####################################################

  list_predictions <- list()
  p = 1


  for (M in uniq_tSk){

    ###################### select dv #############################
    dv_temporal_fit <- subset(df_fit, subset = df_fit$speciesGroup %in% M)

    # BA can not be missing!
    dv_temporal_fit <- dplyr::filter(dv_temporal_fit, !is.na(BA))

    # I change the name of BA into BA_mid, so that the predict function will select the right variable
    dv_temporal_fit <- rename(dv_temporal_fit, "BA_mid" = "BA")

    dv_temporal_predict <- subset(df_predict, subset = df_predict$speciesGroup %in% M)
    # dv_temporal_predict <- dplyr::filter(dv_temporal_predict, !is.na(BA_mid))


    if (height_model == "brnn"){

      if (height_pred_level == 1){

        capture.output(mod1 <- brnn(height ~ BA_mid + plotID, data = dv_temporal_fit, neurons = BRNN_neurons))

      } else {

        capture.output(mod1 <- brnn(height ~ BA_mid, data = dv_temporal_fit, neurons = BRNN_neurons))

      }

    } else {

      mod1 <- lmfor::fithd(dv_temporal_fit$BA_mid,
                    dv_temporal_fit$height,
                    plot = dv_temporal_fit$plotID,
                    modelName = height_model, control=list(maxIter = 1000, msmaxIter = 1000))

    }


    ###################
    # prediction step #
    ###################

    if (height_model == "brnn"){

      # in predict.brnn na.pass is not working, so I do workaround here

      # 1 predictions for BA
      dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))

      dv_temporal_predict_noNA <- dplyr::filter(dv_temporal_predict, !is.na(BA_mid))
      dv_temporal_predict_yesNA <- dplyr::filter(dv_temporal_predict, is.na(BA_mid))
      dv_temporal_predict_noNA$height_mid <- predict(mod1, newdata =dv_temporal_predict_noNA)


      if (sum(is.na(dv_temporal_predict$BA)) > 0){

        dv_temporal_predict_yesNA$height_mid <- NA
        dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dplyr::select(dv_temporal_predict_yesNA, colnames(dv_temporal_predict_noNA)))

      } else {

        dv_temporal_predict <- dv_temporal_predict_noNA

      }

      dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
      dv_temporal_predict$key_temp <- NULL
      dv_temporal_predict$height_new <- NULL
      dv_temporal_predict$p_height_new <- NULL

    } else {

      dv_temporal_predict$height_mid <- predict(mod1, newdata = rename(dv_temporal_predict, "d" = "BA_mid", "plot" = "plotID"), level = height_pred_level, na.action = na.pass)

      dv_temporal_predict$height_new <- NULL
      dv_temporal_predict$p_height_new <- NULL

    }

    list_predictions[[p]] <- dv_temporal_predict
    p = p + 1

  }

  DF_predictions_sGroups <- do.call(rbind, list_predictions)
  DF_predictions_sGroups <-dplyr::filter(DF_predictions_sGroups, !(key %in% DF_predictions_species$key))

  # Now add height data to our data (merge)
  data_height_predictions <- rbind(
    DF_predictions_sGroups,
    DF_predictions_species)

  # correction: Negative height predictions are set to 1 meter (could also be the minimum measured)
  data_height_predictions <- mutate(data_height_predictions, height_mid = ifelse(height_mid < 0, 1, height_mid))

  data_height_predictions$key <- NULL

  return(data_height_predictions)

}

