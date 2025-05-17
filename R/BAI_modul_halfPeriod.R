#' BAI_prediction_halfPeriod
#'
#' BAI model for the MLFS to estimate BAI for half period
#'
#' @return a data frame with calculated basal area increments in the middle of a simulation step
#'
#' @keywords internal

BAI_prediction_halfPeriod <- function(df_fit, df_predict,
                           species_n_threshold = 100,
                           site_vars, include_climate,
                           rf_mtry = NULL,
                           measurement_thresholds = NULL,
                           area_correction = NULL
                           ){

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
  area_factor <- NULL
  weight_mid <- NULL
  ranger <- NULL
  DBH_threshold <- NULL

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
data_below_threshold <- NULL

# Here we define species that will be predicted as individual species
data_above_threshold <- droplevels(df_fit[df_fit$species %in% names(BAI_data)[BAI_data >= species_n_threshold],,drop=FALSE])
unique_dv <- unique(data_above_threshold$species)
data_above_threshold <- NULL

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

    rf_mod <- ranger(formula, data = dv_temporal_fit)

  } else {

    rf_mod <- ranger(formula, data = dv_temporal_fit, mtry = rf_mtry)

  }

  # Do the same for initial data
  dv_temporal_predict <- subset(df_predict, subset = df_predict$species %in% M)

  #### Here we create a work around so we can also apply ranger on data with missing values ###

  dv_temporal_predict_A <- dplyr::filter(dv_temporal_predict, !is.na(BA))
  dv_temporal_predict_B <- dplyr::filter(dv_temporal_predict, is.na(BA))

  temp_predictions<- predict(rf_mod, data = dv_temporal_predict_A)
  dv_temporal_predict_A$BAI_new <- temp_predictions$predictions

  # predicted BAI can't be less than 0
  dv_temporal_predict_A$BAI_new <- ifelse(dv_temporal_predict_A$BAI_new < 0, 0, dv_temporal_predict_A$BAI_new)

  if (nrow(dv_temporal_predict_B) > 0){

    dv_temporal_predict_B$BAI_new <- NA
    dv_temporal_predict <- rbind(dv_temporal_predict_A, dv_temporal_predict_B)

  } else {

    dv_temporal_predict <- dv_temporal_predict_A

  }

  # dv_temporal_predict <- data.table(dv_temporal_predict)
  # class(dv_temporal_predict)
  # dv_temporal_predict[, ':='(BA_mid = BA + BAI_new / 2,
  #                           BAI_mid = BAI_new / 2,
  #                           weight_mid = NA)]

  # dv_temporal_predict[, c("BAI_new"):=NULL]

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
                                       weight_mid = NA)

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

    rf_mod <- ranger(formula, data = dv_temporal_fit)

  } else {

    rf_mod <- ranger(formula, data = dv_temporal_fit, mtry = rf_mtry)
  }

  dv_temporal_predict <- subset(df_predict, subset = df_predict$speciesGroup %in% M)

  #### Here we create a work around so we can also apply ranger on data with missing values ###

  dv_temporal_predict_A <- dplyr::filter(dv_temporal_predict, !is.na(BA))
  dv_temporal_predict_B <- dplyr::filter(dv_temporal_predict, is.na(BA))

  temp_predictions<- predict(rf_mod, data = dv_temporal_predict_A)
  dv_temporal_predict_A$BAI_new <- temp_predictions$predictions
  # predicted BAI can't be less than 0
  dv_temporal_predict_A$BAI_new <- ifelse(dv_temporal_predict_A$BAI_new < 0, 0, dv_temporal_predict_A$BAI_new)

  if (nrow(dv_temporal_predict_B) > 0){

    dv_temporal_predict_B$BAI_new <- NA
    dv_temporal_predict <- rbind(dv_temporal_predict_A, dv_temporal_predict_B)

  } else {

    dv_temporal_predict <- dv_temporal_predict_A

  }

  # temp_predictions<- predict(rf_mod, data = dv_temporal_predict)
  # dv_temporal_predict$BAI_new <- temp_predictions$predictions
  # dv_temporal_predict$BAI_new <- predict(rf_mod, dv_temporal_predict)
  # predicted BAI can't be less than 0
  # dv_temporal_predict$BAI_new <- ifelse(dv_temporal_predict$BAI_new < 0, 0, dv_temporal_predict$BAI_new)

  dv_temporal_predict <- mutate(dv_temporal_predict,
                                # p_BA = BA, p_volume = volume,
                                BA_mid = BA + BAI_new / 2,
                                #year = year,
                                #height = NA, crownHeight = NA,
                                #stand_BA = NA, stand_n = NA, BAL = NA,
                                BAI_mid = BAI_new / 2,
                                BAI_new = NULL, BA_new = NULL,
                                weight_mid = NA) %>%
    filter(!(species %in% unique_dv))

  list_predictions[[p]] <- dv_temporal_predict
  p = p + 1

}

DF_predictions_sGroups <- do.call(rbind, list_predictions)

# DF_predictions_sGroups <- dplyr::filter(DF_predictions_sGroups, !(species %in% unique_dv))

DF_predictions <- rbind(DF_predictions_species, DF_predictions_sGroups)

# Assign the correct weight to

# In case area correction factors are provided we use them to correct plot weights
measurement_thresholds$BA_threshold <- ((measurement_thresholds$DBH_threshold/2)^2 * pi)/10000

#if (!is.null(area_correction)){
#
#  DF_predictions <- DF_predictions %>% mutate(weight_mid = ifelse(BA_mid >= max(measurement_thresholds$BA_threshold),
#                                          measurement_thresholds[, "weight"][which.max(measurement_thresholds$BA_threshold)],
#                                          measurement_thresholds[, "weight"][which.min(measurement_thresholds$BA_threshold)]),
#
#                          DBH_threshold = ifelse(BA_mid >= max(measurement_thresholds$BA_threshold),
#                                                 measurement_thresholds[, "DBH_threshold"][which.max(measurement_thresholds$BA_threshold)],
#                                                 measurement_thresholds[, "DBH_threshold"][which.min(measurement_thresholds$BA_threshold)]))
#
#  DF_predictions <- merge(DF_predictions, area_correction, by = c("plotID", "DBH_threshold"), all.x = TRUE)
#
#  DF_predictions <- dplyr::mutate(DF_predictions, area_factor = ifelse(is.na(area_factor), 1, area_factor),
#                        weight_mid = weight_mid*area_factor, area_factor = NULL, DBH_threshold = NULL)
#
#} else {
#
#  DF_predictions <- DF_predictions %>% mutate(weight_mid = ifelse(BA_mid >= max(measurement_thresholds$BA_threshold),
#                                                                  measurement_thresholds[, "weight"][which.max(measurement_thresholds$BA_threshold)],
#                                                                  measurement_thresholds[, "weight"][which.min(measurement_thresholds$BA_threshold)]))
#}
#
#return(DF_predictions)
#


if (!is.null(area_correction)) {

  # Precompute threshold rows
  max_row <- which.max(measurement_thresholds$BA_threshold)
  min_row <- which.min(measurement_thresholds$BA_threshold)

  high_w    <- measurement_thresholds$weight[max_row]
  low_w     <- measurement_thresholds$weight[min_row]
  high_dbh  <- measurement_thresholds$DBH_threshold[max_row]
  low_dbh   <- measurement_thresholds$DBH_threshold[min_row]

  DF_predictions <- DF_predictions %>%
    mutate(
      weight_mid    = if_else(BA_mid >= measurement_thresholds$BA_threshold[max_row], high_w,   low_w),
      DBH_threshold = if_else(BA_mid >= measurement_thresholds$BA_threshold[max_row], high_dbh, low_dbh)
    )

  DF_predictions <- merge(
    DF_predictions,
    area_correction,
    by = c("plotID", "DBH_threshold"),
    all.x = TRUE
  )

  DF_predictions <- DF_predictions %>%
    mutate(
      area_factor       = if_else(is.na(area_factor), 1, area_factor),
      weight_mid        = weight_mid * area_factor
    ) %>%
    select(-area_factor, -DBH_threshold)

} else {

  # Precompute for simpler branch
  max_row <- which.max(measurement_thresholds$BA_threshold)
  min_row <- which.min(measurement_thresholds$BA_threshold)

  high_w <- measurement_thresholds$weight[max_row]
  low_w  <- measurement_thresholds$weight[min_row]

  DF_predictions <- DF_predictions %>%
    mutate(
      weight_mid = if_else(BA_mid >= measurement_thresholds$BA_threshold[max_row], high_w, low_w)
    )
}

return(DF_predictions)

}

