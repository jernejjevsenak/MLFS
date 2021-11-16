#' BAI_prediction_halfPeriod
#'
#' BAI model for the MLFS to estimate BAI for half period
#' @keywords internal

BAI_prediction_halfPeriod <- function(df_fit, df_predict,
                           species_n_threshold = 100,
                           site_vars, include_climate,
                           rf_mtry = NULL){

  # Define global variables
  code <- NULL
  BA <- NULL
  BA_mid <- NULL
  volume<- NULL
  BAI<- NULL
  BAI_new<- NULL
  year<- NULL
  species<- NULL
  plotID<- NULL
  treeID<- NULL
  speciesGroup<- NULL
  BAI_pred<- NULL
  height <- NULL
  crownHeight <- NULL

#####################################
# 2 Predict BAI for next NFI period #
#####################################

if (include_climate == TRUE){

  site_vars <- c(site_vars, "p_sum", "t_avg")

}

formula <- as.formula(paste0("BAI ~ BA + BAL + stand_BA + stand_n + species +", paste(all_of(site_vars), collapse = "+")))

df_fit <-  dplyr::filter(df_fit, !is.na(BAI))
BAI_data <-  table(df_fit$species)

# Here we define instances that will be predicted using group variable
data_below_threshold <- droplevels(df_fit[df_fit$species %in% names(BAI_data)[BAI_data < species_n_threshold],,drop=FALSE])
uniq_tSk <- unique(data_below_threshold$speciesGroup)

# Here we define species that will be predicted as individual species
data_above_threshold <- droplevels(df_fit[df_fit$species %in% names(BAI_data)[BAI_data >= species_n_threshold],,drop=FALSE])
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

  ####################
  # Prediction phase #
  ####################

  if (is.null(rf_mtry)){

    rf_mod <- randomForest(formula, data = dv_temporal_fit)

  } else {

    rf_mod <- randomForest(formula, data = dv_temporal_fit, mtry = rf_mtry)

  }

  # Do the same for initial data
  dv_temporal_predict <- subset(df_predict, subset = df_predict$species %in% M)
  dv_temporal_predict$BAI_new <- predict(rf_mod, dv_temporal_predict)

  # predicted BAI can't be less than 0
  dv_temporal_predict$BAI_new <- ifelse(dv_temporal_predict$BAI_new < 0, 0, dv_temporal_predict$BAI_new)

  dv_temporal_predict <- mutate(dv_temporal_predict,
                                       #p_BA = BA,
                                       #p_volume = volume,
                                       #p_height = height,
                                       #p_crownHeight = crownHeight,
                                       BA_mid = BA + BAI_new / 2,
                                       # year = year,
                                       # height = NA, crownHeight = NA,
                                       # stand_BA = NA, stand_n = NA, BAL = NA,
                                       BAI_mid = BAI_new / 2,
                                       BAI_new = NULL, BA_new = NULL,
                                       weight_mid = ifelse(BA_mid > 0.07068583, 16.67, 50)
                                )

  list_predictions[[p]] <- dv_temporal_predict
  p = p + 1

}

DF_predictions_species <- do.call(rbind, list_predictions)

##################################################
# Now we repeat the process for minor tree species

list_predictions <- list()
p = 1

for (M in uniq_tSk){

  dv_temporal_fit <- subset(df_fit, subset = df_fit$speciesGroup %in% M)


  ####################
  # Prediction phase #
  ####################

  if (is.null(rf_mtry)){

    rf_mod <- randomForest(formula, data = dv_temporal_fit)

  } else {

    rf_mod <- randomForest(formula, data = dv_temporal_fit, mtry = rf_mtry)

  }

  dv_temporal_predict <- subset(df_predict, subset = df_predict$speciesGroup %in% M)

  dv_temporal_predict$BAI_new <- predict(rf_mod, dv_temporal_predict)

  # predicted BAI can't be less than 0
  dv_temporal_predict$BAI_new <- ifelse(dv_temporal_predict$BAI_new < 0, 0, dv_temporal_predict$BAI_new)

  dv_temporal_predict <- mutate(dv_temporal_predict,
                                # p_BA = BA, p_volume = volume,
                                BA_mid = BA + BAI_new / 2,
                                #year = year,
                                #height = NA, crownHeight = NA,
                                #stand_BA = NA, stand_n = NA, BAL = NA,
                                BAI_mid = BAI_new / 2,
                                BAI_new = NULL, BA_new = NULL,
                                weight_mid = ifelse(BA_mid > 0.07068583, 16.67, 50))

  list_predictions[[p]] <- dv_temporal_predict
  p = p + 1

}


DF_predictions_sGroups <- do.call(rbind, list_predictions)
DF_predictions_sGroups <- dplyr::filter(DF_predictions_sGroups, !(species %in% unique_dv))

DF_predictions <- rbind(DF_predictions_species, DF_predictions_sGroups)

return(DF_predictions)

}

