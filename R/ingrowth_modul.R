#' predict_ingrowth
#'
#' Ingrowth model
#'
#' @keywords internal
#'
predict_ingrowth <- function(df_fit, df_predict, site_vars = site_vars,
                             eval_model_ingrowth = TRUE, k = 10, blocked_cv = TRUE,
                             ingrowth_model = "glm", rf_mtry = NULL){


  # Define Global variables
  species <- NULL
  speciesGroup <- NULL
  year <- NULL
  plotID <- NULL
  stand_BA <- NULL
  stand_n <- NULL
  BAL <- NULL
  ingrowth_small <- NULL
  ingrowth_big <- NULL
  ingrowth_small_pred <- NULL
  ingrowth_big_pred <- NULL
  ing_small <- NULL
  ing_big <- NULL
  DBH <- NULL
  protected <- NULL

  ############
  # Ingrowth #
  ############

  df_before <- df_predict

  sp_group_data <- select(df_before, species, speciesGroup) %>% distinct()

  df_predict <- select(df_predict, year, plotID, stand_BA, stand_n, BAL, all_of(site_vars)) %>%
    group_by(plotID) %>% summarise_all(.funs = mean, na.rm = TRUE)

  # Machine Learning Method for Count Data?
  formula_ing_small = as.formula(paste0("ingrowth_small ~ stand_BA + stand_n + ", paste(site_vars, collapse = "+")))
  formula_ing_big = as.formula(paste0("ingrowth_big ~ stand_BA + stand_n + ", paste(site_vars, collapse = "+")))

  ####################
  # Evaluation phase #
  ####################

  if (eval_model_ingrowth == TRUE){

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

      if (ingrowth_model == "glm"){

        pois_mod_ing_small <- glm(formula_ing_small, data = train, family = "poisson")
        pois_mod_ing_big <- glm(formula_ing_big, data = train, family = "poisson")

        test$ingrowth_small_pred <- round(predict(pois_mod_ing_small, test, type = "response"), 0)
        test$ingrowth_big_pred <- round(predict(pois_mod_ing_big, test, type = "response"), 0)

      } else if (ingrowth_model == "rf"){

        if (is.null(rf_mtry)){

          pois_mod_ing_small <- randomForest(formula_ing_small, data = train)
          pois_mod_ing_big <- randomForest(formula_ing_big, data = train)

        } else {

          pois_mod_ing_small <- randomForest(formula_ing_small, data = train, mtry = rf_mtry)
          pois_mod_ing_big <- randomForest(formula_ing_big, data = train, mtry = rf_mtry)

        }

        test$ingrowth_small_pred <- round(predict(pois_mod_ing_small, test), 0)
        test$ingrowth_big_pred <- round(predict(pois_mod_ing_big, test), 0)

        # Just in case of negative predictions (highly unlikely, but still possible)
        test$ingrowth_small_pred <- ifelse(test$ingrowth_small_pred < 0, 0, test$ingrowth_small_pred)
        test$ingrowth_big_pred <- ifelse(test$ingrowth_big_pred < 0, 0, test$ingrowth_big_pred)

      } else {

        stop("ingrowth_model should be 'glm' or 'rf'" )

      }

      eval_list[[m]] <- test

    }

    df_eval_ingrowth <- do.call(rbind, eval_list)

    df_eval_ingrowth <- select(df_eval_ingrowth, plotID, year,
                               ingrowth_small, ingrowth_big,
                               ingrowth_small_pred, ingrowth_big_pred)
  } else {

    df_eval_ingrowth <- "the argument set_eval_ingrowth is set to FALSE"

  }

  ####################
  # Prediction phase #
  ####################

  if (ingrowth_model == "glm"){

    pois_mod_ing_small <- glm(formula_ing_small, data = df_fit, family = "poisson")
    pois_mod_ing_big <- glm(formula_ing_big, data = df_fit, family = "poisson")

    df_predict$ing_small <- round(predict(pois_mod_ing_small, df_predict, type = "response"), 0)
    df_predict$ing_big <- round(predict(pois_mod_ing_big, df_predict, type = "response"), 0)

  } else if(ingrowth_model == "rf") {

    if (is.null(rf_mtry)){

      pois_mod_ing_small <- randomForest(formula_ing_small, data = df_fit)
      pois_mod_ing_big <- randomForest(formula_ing_big, data = df_fit)

    } else {

      pois_mod_ing_small <- randomForest(formula_ing_small, data = df_fit, mtry = rf_mtry)
      pois_mod_ing_big <- randomForest(formula_ing_big, data = df_fit, mtry = rf_mtry)

    }

    df_predict$ing_small <- round(predict(pois_mod_ing_small, df_predict, type = "response"), 0)
    df_predict$ing_big <- round(predict(pois_mod_ing_big, df_predict, type = "response"), 0)

    # Just in case of negative predictions (highly unlikely, but still possible)
    df_predict$ing_small <- ifelse(df_predict$ing_small < 0, 0, df_predict$ing_small)
    df_predict$ing_big <- ifelse(df_predict$ing_big < 0, 0, df_predict$ing_big)

  } else {

    stop("ingrowth_model should be 'glm' or 'rf'" )

  }

  # Rezultati so OK

  new_trees_small <- dplyr::select(df_predict, plotID, year, all_of(site_vars), ing_small) %>%
    filter(ing_small > 0)

  list_new_trees_small <- list()
  b = 1

  for (i in 1:nrow(new_trees_small)){

    temp <- new_trees_small[i,]

    for (j in 1:temp$ing_small){

      list_new_trees_small[[b]] <- temp
      b = b +1

    }
  }

  new_trees_small <- do.call(rbind, list_new_trees_small)
  new_trees_small$ing_small <- NULL

  x <- runif(nrow(new_trees_small), 10, 15)
  x <- sort(x, decreasing = F)
  probs <- seq(0.5,2, length.out = nrow(new_trees_small))
  x <- x * probs
  x <- x * probs
  x <- x * probs

  x <- scales::rescale(x, to = c(10,15))
  x <- sample(x)

  new_trees_small$DBH <- x
  new_trees_small$code <- 3

  ##############################################################
  # new trees big

  new_trees_big <- dplyr::select(df_predict, plotID, year, all_of(site_vars), ing_big) %>%
    filter(ing_big > 0)

  list_new_trees_big <- list()
  b = 1

  for (i in 1:nrow(new_trees_big)){

    temp <- new_trees_big[i,]

    for (j in 1:temp$ing_big){

      list_new_trees_big[[b]] <- temp
      b = b +1

    }
  }

  new_trees_big <- do.call(rbind, list_new_trees_big)
  new_trees_big$ing_big <- NULL

  x <- runif(nrow(new_trees_big), 30, 40)
  x <- sort(x, decreasing = F)
  probs <- seq(0.5,2, length.out = nrow(new_trees_big))
  x <- x * probs
  x <- x * probs
  x <- x * probs
  x <- scales::rescale(x, to = c(30,40))
  x <- sample(x)

  new_trees_big$DBH <- x
  new_trees_big$code <- 15

  new_trees <- rbind(new_trees_big, new_trees_small)
  new_treeIDs <- sample(seq(max(df_before$treeID) + 1 , length.out = nrow(new_trees)))

  new_trees$treeID <- new_treeIDs

  ###############################
  # treba bo dolo?iti ?e species

  new_trees$species <- NA
  new_trees <- data.frame(new_trees)

  for (o in 1:nrow(new_trees)){

    temp_plID <- new_trees[o, "plotID"]

    species_vector <- as.character(df_before[df_before$plotID == as.character(temp_plID), "species"])
    new_trees[o, "species"] <-  sample(species_vector)[1]

}

  new_trees <- merge(new_trees, sp_group_data, by = "species")
  new_trees <- mutate(new_trees, BA = ((DBH/2)^2 * pi)/10000,
         weight = ifelse(DBH > 30, 16.67, 50))

  new_trees$BAL <- NA
  new_trees$stand_n <- NA
  new_trees$stand_BA <- NA
  new_trees$height <- NA
  new_trees$crownHeight <- NA
  new_trees$BAI <- NA
  new_trees$p_BA <- NA
  new_trees$p_volume <- NA
  new_trees$volume <- NA
  new_trees$p_height <- NA
  new_trees$p_crownHeight <- NA
  # new_trees$protected <- NA

  new_trees$t_avg <- NA
  new_trees$p_sum <- NA

  protected_df <- df_before %>% group_by(plotID) %>% summarise(protected = median(protected))

  new_trees <- merge(new_trees, protected_df, by = "plotID", all.x = TRUE)

  new_trees <-  select(new_trees, colnames(df_before))

  both_df <- rbind(df_before, new_trees)
  both_df$BAL <- NA
  both_df$stand_BA <- NA
  both_df$stand_n <- NA

  final_output_list <- list(

    predicted_ingrowth = both_df,
    eval_ingrowth = df_eval_ingrowth

  )

  return(final_output_list)

}
