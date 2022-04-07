#' height_prediction
#'
#' Height model
#'
#' @param df_fit data frame with tree heights and basal areas for individual trees
#' @param df_predict data frame which will be used for predictions
#' @param species_n_threshold a positive integer defining the minimum number of
#' observations required to treat a species as an independent group
#' @param height_model character string defining the model to be used for height
#' prediction. If 'brnn', then ANN method with Bayesian Regularization is applied.
#' In addition, all 2- and 3- parametric H-D models from lmfor R package are
#' available.
#' @param height_pred_level integer with value 0 or 1 defining the level of
#' prediction for height-diameter (H-D) models. The value 1 defines a plot-level
#' prediction, while the value 0 defines regional-level predictions. Default is
#' 0. If using 1, make sure to have representative plot-level data for each
#' species.
#' @param BRNN_neurons positive integer defining the number of neurons to bo
#' used in the brnn method.
#' @param eval_model_height logical, should the height model be evaluated and
#' returned as the output
#' @param blocked_cv logical, should the blocked cross-validation be used in the
#' evaluation phase?
#' @param k the number of folds to be used in the k fold cross-validation
#'
#' @examples
#' library(MLFS)
#' data(data_tree_heights)
#' data(data_v2)
#'
#' # A) Example with the BRNN method
#' h_predictions <- height_prediction(df_fit = data_tree_heights,
#'                                    df_predict = data_v2,
#'                                    species_n_threshold = 100,
#'                                    height_pred_level = 0,
#'                                    height_model = "brnn",
#'                                    BRNN_neurons = 3,
#'                                    eval_model_height = FALSE,
#'                                    blocked_cv = TRUE, k = 10
#'                                    )
#'
#' predicted_df <- h_predictions$data_height_predictions # df with imputed heights
#' evaluation_df <- h_predictions$data_height_eval # df with evaluation results
#'
#' \dontrun{
#' # B) Example with lmfor
#' # Note: lmfor is currently removed from CRAN so it is also not available in
#' # MLFS
#'
#' library(lmfor)
#' h_predictions <- height_prediction(df_fit = data_tree_heights,
#'                                    df_predict = data_v2,
#'                                    species_n_threshold = 100,
#'                                    height_pred_level = 0,
#'                                    height_model = "naslund",
#'                                    eval_model_height = FALSE,
#'                                    blocked_cv = FALSE, k = 10
#'                                    )
#' }

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

  height_model <- tolower(height_model)

  if (height_model != "brnn"){

    stop("For predicting tree heights, only brnn method is currently available")

  }


  df_fit <- mutate(df_fit, key = paste0(plotID, "_", treeID, "_", year))
  df_predict <- mutate(df_predict, key = paste0(plotID, "_", treeID, "_", year))

  ###################
  # Height prediction
  ###################

  # Remove trees withour BA
  # no_BA <- dplyr::filter(df_predict, is.na(BA))

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

  eval_list <- list()
  m_holder = 1

  for (M in unique_dv){

    ###################### izbira dv ########################################
    dv_temporal_fit <- subset(df_fit, subset = df_fit$species %in% M)
    # BA can not be missing!
    dv_temporal_fit <- dplyr::filter(dv_temporal_fit, !is.na(BA))

    dv_temporal_predict <- dplyr::filter(df_predict, species %in% M)

    if (height_model == "brnn"){

      if (height_pred_level == 1){

        capture.output(mod_species <- brnn(height ~ BA + plotID, data = dv_temporal_fit, neurons = BRNN_neurons))

      } else {

        capture.output(mod_species <- brnn(height ~ BA, data = dv_temporal_fit, neurons = BRNN_neurons))

      }

    } else {

      # mod_species <- lmfor::fithd(dv_temporal_fit$BA,
      #              dv_temporal_fit$height,
      #              plot = dv_temporal_fit$plotID,
      #              modelName = height_model, control=list(maxIter = 1000, msmaxIter = 1000))

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

      eval_step$h_pred <- predict(mod_species, newdata = rename(dv_temporal_fit, "d" = "BA", "plot" = "plotID"), level = height_pred_level, na.action = na.pass)

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

      dv_temporal_predict_noNA <- dplyr::filter(dv_temporal_predict, !is.na(BA))
      dv_temporal_predict_yesNA <- dplyr::filter(dv_temporal_predict, is.na(BA))
      dv_temporal_predict_noNA$height_new <- predict(mod_species, newdata =dv_temporal_predict_noNA)


      if (sum(is.na(dv_temporal_predict$BA)) > 0){

        dv_temporal_predict_yesNA$height_new <- NA
        dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dplyr::select(dv_temporal_predict_yesNA, colnames(dv_temporal_predict_noNA)))

      } else {

        dv_temporal_predict <- dv_temporal_predict_noNA

      }

      dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
      dv_temporal_predict$height <- ifelse(is.na(dv_temporal_predict$height), dv_temporal_predict$height_new, dv_temporal_predict$height)
      dv_temporal_predict$key_temp <- NULL
      dv_temporal_predict$height_new <- NULL
      dv_temporal_predict$p_height_new <- NULL

      # 2 predictions for p_BA
      dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))
      dv_temporal_predict_noNA <- dplyr::filter(dv_temporal_predict, !is.na(p_BA))
      dv_temporal_predict_yesNA <- dplyr::filter(dv_temporal_predict, is.na(p_BA))

      if (nrow(dv_temporal_predict_yesNA) > 0){
        dv_temporal_predict_yesNA$p_height_new <- NA
      }


      dv_temporal_predict_noNA$p_height_new <- predict(mod_species,
                                                   newdata = rename(dv_temporal_predict_noNA,
                                                                          "BA_temp" = "BA",
                                                                          "BA" = "p_BA"))

      if (nrow(dv_temporal_predict_yesNA) > 0){
      dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dplyr::select(dv_temporal_predict_yesNA, colnames(dv_temporal_predict_noNA)))
      } else {
        dv_temporal_predict <- dv_temporal_predict_noNA
      }

      dv_temporal_predict$p_height <- ifelse(is.na(dv_temporal_predict$p_height), dv_temporal_predict$p_height_new, dv_temporal_predict$p_height)

      dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
      dv_temporal_predict$key_temp <- NULL
      dv_temporal_predict$height_new <- NULL
      dv_temporal_predict$p_height_new <- NULL

    } else {

      dv_temporal_predict$height_new <- predict(mod_species, newdata = rename(dv_temporal_predict, "d" = "BA", "plot" = "plotID"), level = height_pred_level, na.action = na.pass)
      dv_temporal_predict$p_height_new <- predict(mod_species, newdata = rename(dv_temporal_predict,
                                                                       "d" = "p_BA", "plot" = "plotID"), level = height_pred_level, na.action = na.pass)

      dv_temporal_predict$height <- ifelse(is.na(dv_temporal_predict$height), dv_temporal_predict$height_new, dv_temporal_predict$height)
      dv_temporal_predict$p_height <- ifelse(is.na(dv_temporal_predict$p_height), dv_temporal_predict$p_height_new, dv_temporal_predict$p_height)

      dv_temporal_predict$height_new <- NULL
      dv_temporal_predict$p_height_new <- NULL

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
    dv_temporal_fit <- dplyr::filter(dv_temporal_fit, !is.na(BA))

    dv_temporal_predict <- subset(df_predict, subset = df_predict$speciesGroup %in% M)
    # dv_temporal_predict <- dplyr::filter(dv_temporal_predict, !is.na(BA))


    if (height_model == "brnn"){

      if (height_pred_level == 1){

        capture.output(mod_speciesGroups <- brnn(height ~ BA + plotID, data = dv_temporal_fit, neurons = BRNN_neurons))

      } else {

        capture.output(mod_speciesGroups <- brnn(height ~ BA, data = dv_temporal_fit, neurons = BRNN_neurons))

      }

    } else {

      # mod_speciesGroups <- lmfor::fithd(dv_temporal_fit$BA,
      #              dv_temporal_fit$height,
      #              plot = dv_temporal_fit$plotID,
      #              modelName = height_model, control=list(maxIter = 1000, msmaxIter = 1000))

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

        eval_step$h_pred <- predict(mod_speciesGroups, newdata = rename(dv_temporal_predict, "d" = "BA", "plot" = "plotID"), level = height_pred_level, na.action = na.pass)

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
      dv_temporal_predict_noNA <- dplyr::filter(dv_temporal_predict, !is.na(BA))
      dv_temporal_predict_yesNA <- dplyr::filter(dv_temporal_predict, is.na(BA))

      dv_temporal_predict_noNA$height_new <- predict(mod_speciesGroups, newdata =dv_temporal_predict_noNA)

    if (sum(is.na(dv_temporal_predict$BA)) > 0){

        dv_temporal_predict_yesNA$height_new <- NA
        dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dplyr::select(dv_temporal_predict_yesNA, colnames(dv_temporal_predict_noNA)))

      } else {

        dv_temporal_predict <- dv_temporal_predict_noNA

      }

      dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
      dv_temporal_predict$height <- ifelse(is.na(dv_temporal_predict$height), dv_temporal_predict$height_new, dv_temporal_predict$height)
      dv_temporal_predict$key_temp <- NULL
      dv_temporal_predict$height_new <- NULL
      dv_temporal_predict$p_height_new <- NULL

      # 2 predictions for p_BA
      dv_temporal_predict$key_temp <- seq(1:nrow(dv_temporal_predict))
      dv_temporal_predict_noNA <- dplyr::filter(dv_temporal_predict, !is.na(p_BA))
      dv_temporal_predict_yesNA <- dplyr::filter(dv_temporal_predict, is.na(p_BA))

      if (nrow(dv_temporal_predict_yesNA) > 0){
        dv_temporal_predict_yesNA$p_height_new <- NA
      }


      dv_temporal_predict_noNA$p_height_new <- predict(mod_speciesGroups,
                                                       newdata = rename(dv_temporal_predict_noNA,
                                                                        "BA_temp" = "BA",
                                                                        "BA" = "p_BA"))

      if (nrow(dv_temporal_predict_yesNA) > 0){
        dv_temporal_predict <- rbind(dv_temporal_predict_noNA, dplyr::select(dv_temporal_predict_yesNA, colnames(dv_temporal_predict_noNA)))
      } else {
        dv_temporal_predict <- dv_temporal_predict_noNA
      }

      dv_temporal_predict$p_height <- ifelse(is.na(dv_temporal_predict$p_height), dv_temporal_predict$p_height_new, dv_temporal_predict$p_height)

      dv_temporal_predict <- arrange(dv_temporal_predict, key_temp)
      dv_temporal_predict$key_temp <- NULL
      dv_temporal_predict$height_new <- NULL
      dv_temporal_predict$p_height_new <- NULL


    } else {

      dv_temporal_predict$height_new <- predict(mod_speciesGroups, newdata = rename(dv_temporal_predict, "d" = "BA", "plot" = "plotID"), level = height_pred_level, na.action = na.pass)
      dv_temporal_predict$p_height_new <- predict(mod_speciesGroups, newdata = rename(dv_temporal_predict,
                                                                       "d" = "p_BA", "plot" = "plotID"), level = height_pred_level, na.action = na.pass)

      dv_temporal_predict$height <- ifelse(is.na(dv_temporal_predict$height), dv_temporal_predict$height_new, dv_temporal_predict$height)
      dv_temporal_predict$p_height <- ifelse(is.na(dv_temporal_predict$p_height), dv_temporal_predict$p_height_new, dv_temporal_predict$p_height)

      dv_temporal_predict$height_new <- NULL
      dv_temporal_predict$p_height_new <- NULL

      }

    list_predictions[[p]] <- dv_temporal_predict
    p = p + 1

  }

  DF_predictions_sGroups <- do.call(rbind, list_predictions)
  DF_predictions_sGroups <-dplyr::filter(DF_predictions_sGroups, !(key %in% DF_predictions_species$key))

  if (eval_model_height == TRUE){

  DF_eval_sGroups <- do.call(rbind, eval_list)
  DF_eval_sGroups <-dplyr::filter(DF_eval_sGroups, !(key %in% DF_eval_species$key))

  data_height_eval <- rbind(
    DF_eval_sGroups,
    DF_eval_species)

  data_height_eval <- dplyr::select(data_height_eval, plotID, treeID, year, speciesGroup, code, species, height, BA, h_pred)

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
    data_height_eval = data_height_eval,
    model_species = mod_species,
    model_speciesGroups = mod_speciesGroups

  )

  return(final_output_list)

}

