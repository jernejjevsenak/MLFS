#' BAI_prediction
#'
#' The Basal Area Increment BAI sub model that is run within the MLFS
#'
#' @param df_fit a data frame with Basal Area Increments (BAI) and all
#' independent variables as specified with the formula
#' @param df_predict data frame which will be used for BAI predictions
#' @param species_n_threshold a positive integer defining the minimum number of
#' observations required to treat a species as an independent group
#' @param site_vars a character vector of variable names which are used as site
#' descriptors
#' @param include_climate logical, should climate variables be included as
#' predictors
#' @param eval_model_BAI logical, should the the BAI model be evaluated and
#' returned as the output
#' @param rf_mtry a number of variables randomly sampled as candidates at
#' each split of a random forest model for predicting basal area increments
#' (BAI). If NULL, default settings are applied.
#' @param k the number of folds to be used in the k fold cross-validation
#' @param blocked_cv logical, should the blocked cross-validation be used in the
#' evaluation phase?
#' @param measurement_thresholds data frame with two variables: 1) DBH_threshold
#' and 2) weight. This information is used to assign the correct weights in BAI
#' and increment sub-model; and to upscale plot-level data to hectares.
#' @param area_correction an optional data frame with three variables: 1) plotID
#' and 2) DBH_threshold and 3) the correction factor to be multiplied by weight
#' for this particular category
#'
#' @examples
#' library(MLFS)
#' data(data_BAI)
#' data(data_v6)
#' data(measurement_thresholds)
#'
#' @return a list with four elements:
#' \enumerate{
#'  \item $predicted_BAI - a data frame with calculated basal area increments (BAI)
#'  \item $eval_BAI - a data frame with predicted and observed basal area increments (BAI), or a character string indicating that BAI model was not evaluated
#'  \item $rf_model_species - the output model for BAI (species level)
#'  \item $rf_model_speciesGroups - the output model for BAI (species group level)
#' }
#'
#' # add BA to measurement thresholds
#' measurement_thresholds$BA_threshold <- ((measurement_thresholds$DBH_threshold/2)^2 * pi)/10000
#'
#' BAI_outputs <- BAI_prediction(df_fit = data_BAI,
#'   df_predict = data_v6,
#'   site_vars = c("slope", "elevation", "northness", "siteIndex"),
#'   rf_mtry = 3,
#'   species_n_threshold = 100,
#'   include_climate = TRUE,
#'   eval_model_BAI = FALSE,
#'   k = 10, blocked_cv = TRUE,
#'   measurement_thresholds = measurement_thresholds)
#'
#' # get the ranger objects
#' BAI_outputs_model_species <- BAI_outputs$rf_model_species
#' BAI_outputs_model_groups <- BAI_outputs$rf_model_speciesGroups
#'

BAI_prediction <- function(df_fit, df_predict,
                           species_n_threshold = 100,
                           site_vars, include_climate,
                           eval_model_BAI = TRUE, rf_mtry = NULL,
                           k = 10, blocked_cv = TRUE,
                           measurement_thresholds = NULL,
                           area_correction = NULL

                           ){

  # Define global variables
  code <- NULL
  BA<- NULL
  volume<- NULL
  BAI<- NULL
  BAI_new<- NULL
  volume_mid<- NULL
  year<- NULL
  species<- NULL
  plotID<- NULL
  treeID<- NULL
  speciesGroup<- NULL
  BAI_pred<- NULL
  height <- NULL
  crownHeight <- NULL
  BA_mid <- NULL
  crownHeight_mid <- NULL
  height_mid <- NULL
  weight <- NULL
  weight_mid <- NULL
  area_factor <- NULL
  ranger <- NULL
  DBH_threshold <- NULL

#####################################
# 2 Predict BAI for next NFI period #
#####################################

dead_trees <- dplyr::filter(df_predict, code %in%  c(1,2)) %>%
                     mutate(p_BA = BA,
                            p_weight = weight,
                            p_height = height,
                            p_crownHeight = crownHeight,
                            p_volume = volume,

                            p_volume_mid = volume_mid,
                            p_BA_mid = BA_mid,
                            p_weight_mid = weight_mid,
                            p_height_mid = height_mid,
                            p_crownHeight_mid = crownHeight_mid,

                            volume_mid = NA,
                            BA_mid = NA,
                            BAI_mid = NA,
                            weight_mid = NA,
                            height_mid = NA,
                            crownHeight_mid = NA,

                            BA = NA,
                            weight = NA,
                            BAI = NA,
                            height = NA,
                            crownHeight = NA)

if (include_climate == TRUE){

  site_vars <- c(site_vars, "p_sum", "t_avg")

}

formula <- as.formula(paste0("BAI ~ BA + BAL + stand_BA + stand_n + species +", paste(all_of(site_vars), collapse = "+")))

df_predict <- dplyr::filter(df_predict, !(code %in% c(1,2)))

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
  unique_species_missing_spG <- unique(df_predict[df_predict$speciesGroup %in% missing_spG,"species"])

  # are these already in unique_dv?
  if (any(!(missing_spG %in% unique_dv))){

    # If not, append it
    uniq_tSk <- c(missing_spG, uniq_tSk)

  }

}

######################################################################
######################################################################

list_predictions <- list()
list_evaluation <- list()
p = 1

eval_list <- list()
m_holder = 1

for (M in unique_dv){

  # select species
  dv_temporal_fit <- subset(df_fit, subset = df_fit$species %in% M)

  ####################
  # Evaluation phase #
  ####################

  if (eval_model_BAI == TRUE){

    if (blocked_cv == FALSE){

      sampled_sequence <- sample(nrow(dv_temporal_fit))

      df_fit_cv <- dv_temporal_fit[sampled_sequence, ]

    } else {

      df_fit_cv <- dv_temporal_fit

    }

    foldi <- seq(1:k)
    folds <- cut(seq(1, nrow(dv_temporal_fit)), breaks = k, labels = FALSE)

    for (m in 1:k){

      testIndexes <- which(folds == m, arr.ind = TRUE)

      test <- df_fit_cv[testIndexes, ]
      train <- df_fit_cv[-testIndexes, ]

      if (is.null(rf_mtry)){

        rf_mod_species <- ranger(formula, data = train)

      } else {

        rf_mod_species <- ranger(formula, data = train, mtry = rf_mtry)

      }


      test$BAI_pred <- predict(rf_mod_species, test)$predictions

      eval_list[[m_holder]] <- test
      m_holder = m_holder + 1
    }

  }

  ####################
  # Prediction phase #
  ####################

  if (is.null(rf_mtry)){

    rf_mod_species <- ranger(formula, data = dv_temporal_fit)

  } else {

    rf_mod_species <- ranger(formula, data = dv_temporal_fit, mtry = rf_mtry)

  }

  # Do the same for initial data
  dv_temporal_predict <- subset(df_predict, subset = df_predict$species %in% M)
  dv_temporal_predict$BAI_new <- predict(rf_mod_species, dv_temporal_predict)$predictions

  # predicted BAI can't be less than 0
  dv_temporal_predict$BAI_new <- ifelse(dv_temporal_predict$BAI_new < 0, 0, dv_temporal_predict$BAI_new)

  dv_temporal_predict <- mutate(dv_temporal_predict,
                                       p_BA = BA,
                                       p_volume = volume,
                                       p_height = height,
                                       p_weight = weight,
                                       p_crownHeight = crownHeight,

                                       p_volume_mid = volume_mid,
                                       p_BA_mid = BA_mid,
                                       p_weight_mid = weight_mid,
                                       p_height_mid = height_mid,
                                       p_crownHeight_mid = crownHeight_mid,

                                       volume_mid = NA,
                                       BA_mid = NA,
                                       BAI_mid = NA,
                                       weight_mid = NA,
                                       height_mid = NA,
                                       crownHeight_mid = NA,

                                       BA = BA + BAI_new, year = year,
                                       height = NA,
                                       crownHeight = NA,
                                       stand_BA = NA, stand_n = NA, BAL = NA, BAI = BAI_new,
                                       BAI_new = NULL, BA_new = NULL,
                                       weight = NA)

  list_predictions[[p]] <- dv_temporal_predict
  p = p + 1

}

DF_evaluation_species <- do.call(rbind, eval_list)

DF_predictions_species <- do.call(rbind, list_predictions)

##################################################
# Now we repeat the process for minor tree species

list_predictions <- list()
p = 1

eval_list <- list()
m_holder = 1

for (M in uniq_tSk){

  dv_temporal_fit <- subset(df_fit, subset = df_fit$speciesGroup %in% M)

  ####################
  # Evaluation phase #
  ####################

  if (eval_model_BAI == TRUE){

    if (blocked_cv == FALSE){

      sampled_sequence <- sample(nrow(dv_temporal_fit))

      df_fit_cv <- dv_temporal_fit[sampled_sequence, ]

    } else {

      df_fit_cv <- dv_temporal_fit

    }

    foldi <- seq(1:k)
    folds <- cut(seq(1, nrow(dv_temporal_fit)), breaks = k, labels = FALSE)

    for (m in 1:k){

      testIndexes <- which(folds == m, arr.ind = TRUE)

      test <- df_fit_cv[testIndexes, ]
      train <- df_fit_cv[-testIndexes, ]

      if (is.null(rf_mtry)){

        rf_mod_speciesGroups <- ranger(formula, data = train)

      } else {

        rf_mod_speciesGroups <- ranger(formula, data = train, mtry = rf_mtry)

      }

      test$BAI_pred <- predict(rf_mod_speciesGroups, test)$predictions

      eval_list[[m_holder]] <- test
      m_holder = m_holder + 1
    }

  }

  ####################
  # Prediction phase #
  ####################

  if (is.null(rf_mtry)){

    rf_mod_speciesGroups <- ranger(formula, data = dv_temporal_fit)

  } else {

    rf_mod_speciesGroups <- ranger(formula, data = dv_temporal_fit, mtry = rf_mtry)

  }

  dv_temporal_predict <- subset(df_predict, subset = df_predict$speciesGroup %in% M)

  dv_temporal_predict$BAI_new <- predict(rf_mod_speciesGroups, dv_temporal_predict)$predictions
  # predicted BAI can't be less than 0
  dv_temporal_predict$BAI_new <- ifelse(dv_temporal_predict$BAI_new < 0, 0, dv_temporal_predict$BAI_new)

  dv_temporal_predict <- mutate(dv_temporal_predict,
                                p_BA = BA,
                                p_weight = weight,
                                p_volume = volume,
                                BA = BA + BAI_new,
                                year = year,

                                p_volume_mid = volume_mid,
                                p_BA_mid = BA_mid,
                                p_weight_mid = weight_mid,
                                p_height_mid = height_mid,
                                p_crownHeight_mid = crownHeight_mid,

                                volume_mid = NA,
                                BA_mid = NA,
                                BAI_mid = NA,
                                weight_mid = NA,
                                height_mid = NA,
                                crownHeight_mid = NA,

                                height = NA, crownHeight = NA,
                                stand_BA = NA, stand_n = NA, BAL = NA, BAI = BAI_new,
                                BAI_new = NULL, BA_new = NULL,
                                weight = NA) %>%
    filter(!(species %in% unique_dv))

  list_predictions[[p]] <- dv_temporal_predict
  p = p + 1

}

if (eval_model_BAI == TRUE){

DF_evaluation_sGroups <- do.call(rbind, eval_list)

# DF_evaluation_sGroups <- dplyr::filter(DF_evaluation_sGroups, !(species %in% unique_dv))

data_eval_BAI <- rbind(DF_evaluation_species, DF_evaluation_sGroups)

data_eval_BAI <- dplyr::select(data_eval_BAI, plotID, treeID, year, speciesGroup, code,
                                            species, BAI, BAI_pred)

} else {

  data_eval_BAI <- "the argument set_eval_BAI is set to FALSE"

}

DF_predictions_sGroups <- do.call(rbind, list_predictions)
DF_predictions_sGroups <- dplyr::filter(DF_predictions_sGroups, !(species %in% unique_dv))

# rbind predictions
dead_trees <- dplyr::select(dead_trees, colnames(DF_predictions_species))
DF_predictions <- rbind(dead_trees, DF_predictions_species, DF_predictions_sGroups)


# Assign the correct weight to

# In case area correction factors are provided we use them to correct plot weights
measurement_thresholds$BA_threshold <- ((measurement_thresholds$DBH_threshold/2)^2 * pi)/10000







#if (!is.null(area_correction)){
#
#  DF_predictions <- DF_predictions %>% mutate(weight = ifelse(BA >= max(measurement_thresholds$BA_threshold),
#                                                                  measurement_thresholds[, "weight"][which.max(measurement_thresholds$BA_threshold)],
#                                                                  measurement_thresholds[, "weight"][which.min(measurement_thresholds$BA_threshold)]),
#
#                                              DBH_threshold = ifelse(BA >= max(measurement_thresholds$BA_threshold),
#                                                                     measurement_thresholds[, "DBH_threshold"][which.max(measurement_thresholds$BA_threshold)],
#                                                                     measurement_thresholds[, "DBH_threshold"][which.min(measurement_thresholds$BA_threshold)]))
#
#  DF_predictions <- merge(DF_predictions, area_correction, by = c("plotID", "DBH_threshold"), all.x = TRUE)
#
#  DF_predictions <- dplyr::mutate(DF_predictions, area_factor = ifelse(is.na(area_factor), 1, area_factor),
#                                  weight = weight*area_factor, area_factor = NULL, DBH_threshold = NULL)
#
#} else {
#
#  DF_predictions <- DF_predictions %>% mutate(weight = ifelse(BA >= max(measurement_thresholds$BA_threshold),
#                                                                  measurement_thresholds[, "weight"][which.max(measurement_thresholds$BA_threshold)],
#                                                                  measurement_thresholds[, "weight"][which.min(measurement_thresholds$BA_threshold)]))
#}

#final_output_list <- list(
#
#  predicted_BAI = DF_predictions,
#  eval_BAI = data_eval_BAI,
#  rf_model_species = rf_mod_species,
#  rf_model_speciesGroups = rf_mod_speciesGroups
#
#)

if (!is.null(area_correction)) {

  # Precompute threshold rows and corresponding values
  max_row <- which.max(measurement_thresholds$BA_threshold)
  min_row <- which.min(measurement_thresholds$BA_threshold)

  high_w    <- measurement_thresholds$weight[max_row]
  low_w     <- measurement_thresholds$weight[min_row]
  high_dbh  <- measurement_thresholds$DBH_threshold[max_row]
  low_dbh   <- measurement_thresholds$DBH_threshold[min_row]

  DF_predictions <- DF_predictions %>%
    mutate(
      weight        = if_else(BA >= measurement_thresholds$BA_threshold[max_row], high_w,   low_w),
      DBH_threshold = if_else(BA >= measurement_thresholds$BA_threshold[max_row], high_dbh, low_dbh)
    )

  DF_predictions <- merge(
    DF_predictions,
    area_correction,
    by = c("plotID", "DBH_threshold"),
    all.x = TRUE
  )

  DF_predictions <- DF_predictions %>%
    mutate(
      area_factor = if_else(is.na(area_factor), 1, area_factor),
      weight      = weight * area_factor
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
      weight = if_else(BA >= measurement_thresholds$BA_threshold[max_row], high_w, low_w)
    )
}

final_output_list <- list(
  predicted_BAI         = DF_predictions,
  eval_BAI              = data_eval_BAI,
  rf_model_species      = rf_mod_species,
  rf_model_speciesGroups = rf_mod_speciesGroups
)

return(final_output_list)

}

