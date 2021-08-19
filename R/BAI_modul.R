#' BAI_prediction
#'
#' BAI model for the MLFS
#' @keywords internal

BAI_prediction <- function(df_fit, df_predict,
                           species_n_threshold = 100,
                           site_vars, include_climate,
                           eval_model_BAI = TRUE, rf_mtry = NULL,
                           k = 10, blocked_cv = TRUE

                           ){

  # Define global variables
  code <- NULL
  BA<- NULL
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

dead_trees <-filter(df_predict, code %in%  c(1,2)) %>%
    mutate(p_BA = BA,
           p_height = height,
           p_volume = volume,
           BA = NA,
           BAI = NA,
           height = NA,
           crownHeight = NA)

if (include_climate == TRUE){

  site_vars <- c(site_vars, "p_sum", "t_avg")

}

formula <- as.formula(paste0("BAI ~ BA + BAL + stand_BA + stand_n + species +", paste(site_vars, collapse = "+")))

df_predict <-filter(df_predict, !(code %in% c(1,2)))

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

        rf_mod <- randomForest(formula, data = train)

      } else {

        rf_mod <- randomForest(formula, data = train, mtry = rf_mtry)

      }

      test$BAI_pred <- predict(rf_mod, test)

      eval_list[[m_holder]] <- test
      m_holder = m_holder + 1
    }

  }

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

  dv_temporal_predict <- mutate(dv_temporal_predict,
                                       p_BA = BA,
                                       p_volume = volume,
                                       p_height = height,
                                       p_crownHeight = crownHeight,
                                       BA = BA + BAI_new, year = year,
                                       height = NA, crownHeight = NA,
                                       stand_BA = NA, stand_n = NA, BAL = NA, BAI = BAI_new,
                                       BAI_new = NULL, BA_new = NULL,
                                       weight = ifelse(BA > 0.07068583, 16.67, 50))

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

        rf_mod <- randomForest(formula, data = train)

      } else {

        rf_mod <- randomForest(formula, data = train, mtry = rf_mtry)

      }

      test$BAI_pred <- predict(rf_mod, test)

      eval_list[[m_holder]] <- test
      m_holder = m_holder + 1
    }

  }

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

  dv_temporal_predict <- mutate(dv_temporal_predict,
                                p_BA = BA, p_volume = volume, BA = BA + BAI_new, year = year,
                            height = NA, crownHeight = NA,
                            stand_BA = NA, stand_n = NA, BAL = NA, BAI = BAI_new,
                            BAI_new = NULL, BA_new = NULL,
                            weight = ifelse(BA > 0.07068583, 16.67, 50))

  list_predictions[[p]] <- dv_temporal_predict
  p = p + 1

}

if (eval_model_BAI == TRUE){

DF_evaluation_sGroups <- do.call(rbind, eval_list)
DF_evaluation_sGroups <- filter(DF_evaluation_sGroups, !(species %in% unique_dv))

data_eval_BAI <- rbind(DF_evaluation_species, DF_evaluation_sGroups)

data_eval_BAI <- select(data_eval_BAI, plotID, treeID, year, speciesGroup, code,
                                            species, BAI, BAI_pred)

} else {

  data_eval_BAI <- "the argument set_eval_BAI is set to FALSE"

}

DF_predictions_sGroups <- do.call(rbind, list_predictions)
DF_predictions_sGroups <- filter(DF_predictions_sGroups, !(species %in% unique_dv))

# rbind predictions
dead_trees <- select(dead_trees, colnames(DF_predictions_species))
DF_predictions <- rbind(dead_trees, DF_predictions_species, DF_predictions_sGroups)

final_output_list <- list(

  predicted_BAI = DF_predictions,
  eval_BAI = data_eval_BAI

)

return(final_output_list)

}

