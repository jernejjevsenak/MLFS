#' crownHeight_prediction_halfPeriod
#'
#' Models to predict crown height (half period)
#'
#' @return a data frame with imputed crown heights in the middle of a simulation step
#'
#' @keywords internal
#'
crownHeight_prediction_halfPeriod <- function(df_fit,  df_predict,
                                      site_vars = site_vars,
                                      species_n_threshold = 100,
                                      crownHeight_model = "lm",
                                      BRNN_neurons = 3){

  # Define global variables
  plotID <- NULL
  treeID <- NULL
  year <- NULL
  BA <- NULL
  crownHeight <- NULL
  crownHeight_mid <- NULL
  height <- NULL
  height_mid <- NULL
  p_height <- NULL
  key_temp <- NULL
  key <- NULL
  speciesGroup <- NULL
  code <- NULL
  species <- NULL
  crownHeight_pred <- NULL

  df_fit <- mutate(df_fit, key = paste0(plotID, "_", treeID, "_", year))
  df_predict <- mutate(df_predict, key = paste0(plotID, "_", treeID, "_", year))

  # We define here strategy: if there are enough measurements, we calculate crownHeights for species
  # If not, we use grouping ID
  crownHeight_data <-  dplyr::filter(df_fit, !is.na(crownHeight))
  crownHeight_data <-  table(crownHeight_data$species)

  # Here we define instances that will be predicted using group variable
  data_below_threshold <- droplevels(df_fit[df_fit$species %in% names(crownHeight_data)[crownHeight_data < species_n_threshold],,drop=FALSE])
  uniq_tSk <- unique(data_below_threshold$speciesGroup)

  # Here we define species that will be predicted as individual species
  data_above_threshold <- droplevels(df_fit[df_fit$species %in% names(crownHeight_data)[crownHeight_data >= species_n_threshold],,drop=FALSE])
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

    # select species
    dv_temporal_fit <- subset(df_fit, subset = df_fit$species %in% M)

    # BA can not be missing!
    dv_temporal_fit <- dplyr::filter(dv_temporal_fit, !is.na(BA))

    # I change the name of BA into BA_mid, so that the predict function will select the right variable
    dv_temporal_fit <- rename(dv_temporal_fit, "BA_mid" = "BA", "height_mid" = "height")

    temp_df <- dplyr::select(dv_temporal_fit, crownHeight, height_mid, all_of(site_vars))

    dv_temporal_predict <- subset(df_predict, subset = df_predict$species %in% M)




    ###################
    # prediction part #
    ###################

    if (crownHeight_model == "lm"){

      mod1 <- lm(crownHeight ~ ., data = temp_df)

    } else if(crownHeight_model == "brnn") {

      capture.output(mod1 <- brnn(crownHeight ~ ., data = temp_df, neurons = BRNN_neurons))

    }

    # predict crownHeight for BA
    dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))
    dv_temporal_predict_noNA <- dplyr::filter(dv_temporal_predict, !is.na(height_mid))
    dv_temporal_predict_yesNA <- dplyr::filter(dv_temporal_predict, is.na(height_mid))
    dv_temporal_predict_noNA$crownHeight_mid <- predict(mod1, newdata = dv_temporal_predict_noNA)


    if (sum(is.na(dv_temporal_predict$BA)) > 0){

      dv_temporal_predict_yesNA$crownHeight_mid <- NA
      dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dplyr::select(dv_temporal_predict_yesNA, colnames(dv_temporal_predict_noNA)))

    } else {

      dv_temporal_predict <- dv_temporal_predict_noNA

    }

    # dv_temporal_predict$crownHeight <- ifelse(is.na(dv_temporal_predict$crownHeight), dv_temporal_predict$crownHeight_new, dv_temporal_predict$crownHeight)

    dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
    dv_temporal_predict$key_temp <- NULL
    dv_temporal_predict$crownHeight_new <- NULL

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

    # select species group
    dv_temporal_fit <- subset(df_fit, subset = df_fit$speciesGroup %in% M)

    # BA can not be missing!
    dv_temporal_fit <- dplyr::filter(dv_temporal_fit, !is.na(BA))

    # I change the name of BA into BA_mid, so that the predict function will select the right variable
    dv_temporal_fit <- rename(dv_temporal_fit, "BA_mid" = "BA", "height_mid" = "height")

    temp_df <- dplyr::select(dv_temporal_fit, crownHeight, height_mid, all_of(site_vars))

    dv_temporal_predict <- subset(df_predict, subset = df_predict$speciesGroup %in% M)

    ###################
    # prediction part #
    ###################

    if (crownHeight_model == "lm"){

      mod1 <- lm(crownHeight ~ ., data = temp_df)

    } else if(crownHeight_model == "brnn") {

      capture.output(mod1 <- brnn(crownHeight ~ ., data = temp_df, neurons = BRNN_neurons))

    }

    # predict crownHeight for BA
    dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))

    dv_temporal_predict_noNA <- dplyr::filter(dv_temporal_predict, !is.na(height_mid))
    dv_temporal_predict_yesNA <- dplyr::filter(dv_temporal_predict, is.na(height_mid))
    dv_temporal_predict_noNA$crownHeight_mid <- predict(mod1, newdata = dv_temporal_predict_noNA)

    if (sum(is.na(dv_temporal_predict$BA)) > 0){

      dv_temporal_predict_yesNA$crownHeight_mid <- NA
      dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dv_temporal_predict_yesNA)

    } else {

      dv_temporal_predict <- dv_temporal_predict_noNA

    }


    # dv_temporal_predict$crownHeight <- ifelse(is.na(dv_temporal_predict$crownHeight), dv_temporal_predict$crownHeight_new, dv_temporal_predict$crownHeight)
    dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
    dv_temporal_predict$key_temp <- NULL
    dv_temporal_predict$crownHeight_new <- NULL

    list_predictions[[p]] <- dv_temporal_predict
    p = p + 1

  }

  DF_predictions_sGroups <- do.call(rbind, list_predictions)
  DF_predictions_sGroups <-dplyr::filter(DF_predictions_sGroups, !(key %in% DF_predictions_species$key))


  data_crownHeight_predictions <- rbind(
    DF_predictions_sGroups,
    DF_predictions_species)

  # correction: Negative crownHeight predictions are set to 1 meter (could also be the minimum measured)
  data_crownHeight_predictions <- mutate(data_crownHeight_predictions, crownHeight_mid = ifelse(crownHeight_mid < 0, 1, crownHeight_mid))

  data_crownHeight_predictions$key <- NULL

  return(data_crownHeight_predictions)

}

