#' predict_ingrowth
#'
#' Ingrowth model for predicting new trees witihn the MLFS
#'
#' @param df_fit a plot-level data with plotID, stand variables and site
#' descriptors, and the two target variables describing the number of ingrowth
#' trees for inner (ingrowth_3) and outer (ingrowth_15) circles
#' @param df_predict data frame which will be used for ingrowth predictions
#' @param site_vars a character vector of variable names which are used as site
#' descriptors
#' @param include_climate logical, should climate variables be included as
#' predictors
#' @param eval_model_ingrowth logical, should the the ingrowth model be
#' evaluated and returned as the output
#' @param rf_mtry a number of variables randomly sampled as candidates at
#' each split of a random forest model for predicting ingrowth. If NULL, default
#' settings are applied.
#' @param k the number of folds to be used in the k fold cross-validation
#' @param blocked_cv logical, should the blocked cross-validation be used in the
#' evaluation phase?
#' @param ingrowth_model model to be used for ingrowth predictions. 'glm' for
#' generalized linear models (Poisson regression), 'ZIF_poiss' for zero inflated
#' Poisson regression and 'rf' for random forest
#' @param ingrowth_table a data frame with 4 variables: (ingrowth) code,
#' DBH_threshold, DBH_max and weight. Ingrowth table is used within the ingrowth
#' sub model to correctly simulate different ingrowth levels and associated
#' upscale weights
#' @param DBH_distribution_parameters A list with deciles of DBH distributions
#' that are used to simulate DBH for new trees, seperately for each ingrowth
#' category
#'
#' @examples
#'
#' library(MLFS)
#'
#' data(data_v6)
#' data(data_ingrowth)
#' data(ingrowth_table)
#' data(ingrowth_parameter_list)
#'
#' ingrowth_outputs <- predict_ingrowth(
#'    df_fit = data_ingrowth,
#'    df_predict = data_v6,
#'    site_vars = c("slope", "elevation", "northness", "siteIndex"),
#'    include_climate = TRUE,
#'    eval_model_ingrowth = TRUE,
#'    rf_mtry = 3,
#'    k = 10, blocked_cv = TRUE,
#'    ingrowth_model = 'rf',
#'    ingrowth_table = ingrowth_table,
#'    DBH_distribution_parameters = ingrowth_parameter_list)
#'

predict_ingrowth <- function(df_fit, df_predict, site_vars = site_vars,
                             include_climate = include_climate,
                             eval_model_ingrowth = TRUE, k = 10, blocked_cv = TRUE,
                             ingrowth_model = "glm", rf_mtry = NULL,
                             ingrowth_table = NULL, DBH_distribution_parameters = NULL){

  # Define Global variables
  species <- NULL
  speciesGroup <- NULL
  year <- NULL
  plotID <- NULL
  stand_BA <- NULL
  stand_n <- NULL
  BAL <- NULL
  ing_small <- NULL
  ing_big <- NULL
  DBH <- NULL
  protected <- NULL
  code <- NULL
  threshold <- NULL
  weight <- NULL

  ############
  # Ingrowth #
  ############

  df_before <- df_predict

  if (include_climate == TRUE){

    site_vars_A <- c(site_vars, "p_sum", "t_avg")

  } else {

    site_vars_A <- site_vars

  }

  sp_group_data <- dplyr::select(df_before, species, speciesGroup) %>% distinct()

  df_predict <- dplyr::select(df_predict, year, plotID, stand_BA, stand_n, BAL, all_of(site_vars_A)) %>%
    group_by(plotID) %>% summarise_all(.funs = mean, na.rm = TRUE)


  # extract the
  ing_codes <- unique(ingrowth_table$code)

  for (i_codes in ing_codes){

    assign(paste0("formula_ing_", i_codes), as.formula(paste(paste0("ingrowth_",i_codes), "~ stand_BA + stand_n + ", paste(site_vars, collapse = "+"))))

  }

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

        for (i_codes in ing_codes){

         assign(paste0("mod_ing_", i_codes), eval(parse(text = paste0("glm(formula_ing_", i_codes, ", data = train, family = 'poisson')"))))
         assign("new_temp_var", eval(parse(text = paste0("round(predict(mod_ing_", i_codes, ", test, type = 'response'), 0)"))))
         test$new_temp_var  <- new_temp_var
         colnames(test)[ncol(test)] <- paste0("ingrowth_",i_codes,"_pred")

         }

      } else if  (ingrowth_model == "ZIF_poiss"){

        for (i_codes in ing_codes){

          assign(paste0("mod_ing_", i_codes), eval(parse(text = paste0("zeroinfl(formula_ing_", i_codes, ", data = train)"))))

          assign("new_temp_var", eval(parse(text = paste0("round(predict(mod_ing_", i_codes, ", test, type = 'response'), 0)"))))
          test$new_temp_var  <- new_temp_var
          colnames(test)[ncol(test)] <- paste0("ingrowth_",i_codes,"_pred")

        }

      } else if (ingrowth_model == "rf"){

        for (i_codes in ing_codes){

        if (is.null(rf_mtry)){

          assign(paste0("mod_ing_", i_codes), eval(parse(text = paste0("ranger(formula_ing_", i_codes, ", data = train)"))))
          assign("new_temp_var", eval(parse(text = paste0("round(predict(mod_ing_", i_codes, ", test, type = 'response')$predictions, 0)"))))
          new_temp_var <- ifelse(new_temp_var < 0, 0, new_temp_var) # In case of negative predictions
          test$new_temp_var  <- new_temp_var
          colnames(test)[ncol(test)] <- paste0("ingrowth_",i_codes,"_pred")

        } else {

          assign(paste0("mod_ing_", i_codes), eval(parse(text = paste0("ranger(formula_ing_", i_codes, ", data = train, mtry = ",rf_mtry,")"))))
          assign("new_temp_var", eval(parse(text = paste0("round(predict(mod_ing_", i_codes, ", test, type = 'response')$predictions, 0)"))))
          new_temp_var <- ifelse(new_temp_var < 0, 0, new_temp_var) # In case of negative predictions
          test$new_temp_var  <- new_temp_var
          colnames(test)[ncol(test)] <- paste0("ingrowth_",i_codes,"_pred")
        }

      }

      } else {

        stop("ingrowth_model should be 'glm', 'ZIF_poiss', or 'rf'" )

      }

      eval_list[[m]] <- test

    }

    df_eval_ingrowth <- do.call(rbind, eval_list)

    c(paste0("ingrowth_", ing_codes), paste0("ingrowth_", ing_codes, "_pred"))

    df_eval_ingrowth <- dplyr::select(df_eval_ingrowth, plotID, year,
                                      c(paste0("ingrowth_", ing_codes),
                                        paste0("ingrowth_", ing_codes, "_pred")))
  } else {

    df_eval_ingrowth <- "the argument set_eval_ingrowth is set to FALSE"

  }

  ####################
  # Prediction phase #
  ####################

  if (ingrowth_model == "glm"){

    for (i_codes in ing_codes){

      assign(paste0("mod_ing_", i_codes), eval(parse(text = paste0("glm(formula_ing_", i_codes, ", data = df_fit, family = 'poisson')"))))

      assign("new_temp_var", eval(parse(text = paste0("round(predict(mod_ing_", i_codes, ", df_predict, type = 'response'), 0)"))))
      df_predict$new_temp_var  <- new_temp_var
      colnames(df_predict)[ncol(df_predict)] <- paste0("ingrowth_",i_codes)

    }


  } else if  (ingrowth_model == "ZIF_poiss"){

    for (i_codes in ing_codes){

      assign(paste0("mod_ing_", i_codes), eval(parse(text = paste0("zeroinfl(formula_ing_", i_codes, ", data = df_fit)"))))
      assign("new_temp_var", eval(parse(text = paste0("round(predict(mod_ing_", i_codes, ", df_predict, type = 'response'), 0)"))))
      df_predict$new_temp_var  <- new_temp_var
      colnames(df_predict)[ncol(df_predict)] <- paste0("ingrowth_",i_codes)

    }

 } else if(ingrowth_model == "rf") {

   for (i_codes in ing_codes){

     if (is.null(rf_mtry)){

       assign(paste0("mod_ing_", i_codes), eval(parse(text = paste0("ranger(formula_ing_", i_codes, ", data = df_fit)"))))
       assign("new_temp_var", eval(parse(text = paste0("round(predict(mod_ing_", i_codes, ", df_predict, type = 'response')$predictions, 0)"))))
       new_temp_var <- ifelse(new_temp_var < 0, 0, new_temp_var) # In case of negative predictions
       df_predict$new_temp_var  <- new_temp_var
       colnames(df_predict)[ncol(df_predict)] <- paste0("ingrowth_",i_codes)

     } else {

       assign(paste0("mod_ing_", i_codes), eval(parse(text = paste0("ranger(formula_ing_", i_codes, ", data = df_fit, mtry = ",rf_mtry,")"))))
       assign("new_temp_var", eval(parse(text = paste0("round(predict(mod_ing_", i_codes, ", df_predict, type = 'response')$predictions, 0)"))))
       new_temp_var <- ifelse(new_temp_var < 0, 0, new_temp_var) # In case of negative predictions
       df_predict$new_temp_var  <- new_temp_var
       colnames(df_predict)[ncol(df_predict)] <- paste0("ingrowth_",i_codes)
     }

   }

  } else {

    stop("ingrowth_model should be 'glm', 'ZIF_poiss', or 'rf'" )

  }

  # Creating the DBH distributions of new trees
  new_trees_list_bind <- list()
  b = 1

  for (i_codes in ing_codes){

    assign("new_trees_temp", eval(parse(text = paste0("dplyr::select(df_predict, plotID, year, all_of(site_vars_A),ingrowth_",i_codes,") %>% dplyr::filter(ingrowth_", i_codes," > 0)"))))

    list_new_trees <- list()
    b = 1

    for (i in 1:nrow(new_trees_temp)){

      temp <- new_trees_temp[i,]

      row_n <- as.numeric(get("temp")[paste0("ingrowth_", i_codes)])
      row_n <- ifelse(row_n > 150, 150, row_n)

      for (j in 1:row_n){

        list_new_trees[[b]] <- temp
        b = b +1

      }
    }

    assign("new_trees_temp", do.call(rbind, list_new_trees))

    assign("new_trees_temp", eval(parse(text = paste0("dplyr::select(new_trees_temp, -ingrowth_",i_codes, ")"))))

    temp_parameters <- unname(DBH_distribution_parameters[which(names(DBH_distribution_parameters) == i_codes)] )[[1]]

    nn <- nrow(new_trees_temp) /  (length(temp_parameters) - 1) + 1

    q_list <- list()

    for (i in 1:(length(temp_parameters) - 1)){

      temp_sim <- runif(nn,temp_parameters[i],temp_parameters[i + 1])
      q_list[[i]] <- temp_sim

    }

    sim_distribution <- c(do.call(rbind, q_list))

    # due to rounding, the numbers usually won't match: apply correction
    if (nrow(new_trees_temp) != length(sim_distribution)){

      diff_n <- length(sim_distribution) - nrow(new_trees_temp)

      sim_distribution <- sim_distribution[-sample(length(sim_distribution), diff_n)]
    }

    new_trees_temp$DBH <- sim_distribution
    new_trees_temp$code <- i_codes

    new_trees_list_bind[[b]] <- new_trees_temp
    b = b + 1
  }

  new_trees <- do.call(rbind, new_trees_list_bind)
  new_treeIDs <- sample(seq(max(df_before$treeID) + 1 , length.out = nrow(new_trees)))

  new_trees$treeID <- new_treeIDs

  ###############################
  # Assigning a species code to new trees

  new_trees$species <- NA
  new_trees <- data.frame(new_trees)

  for (o in 1:nrow(new_trees)){

    temp_plID <- new_trees[o, "plotID"]

    species_vector <- as.character(df_before[df_before$plotID == as.character(temp_plID), "species"])
    new_trees[o, "species"] <-  sample(species_vector)[1]

}

  new_trees <- merge(new_trees, sp_group_data, by = "species")
  new_trees <- mutate(new_trees, BA = ((DBH/2)^2 * pi)/10000)

  new_trees <- merge(new_trees, dplyr::select(ingrowth_table, code, weight), by = "code", all.x = TRUE)

  new_trees$BAL <- NA
  new_trees$stand_n <- NA
  new_trees$stand_BA <- NA
  new_trees$height <- NA
  new_trees$crownHeight <- NA
  new_trees$BAI <- NA
  new_trees$p_BA <- NA
  new_trees$p_weight <- NA
  new_trees$p_volume <- NA
  new_trees$p_volume_mid <- NA
  new_trees$volume <- NA
  new_trees$p_height <- NA
  new_trees$p_crownHeight <- NA
  new_trees$p_height <- NA

  new_trees$p_BA_mid <- NA
  new_trees$p_weight_mid <- NA
  new_trees$p_height_mid <- NA
  new_trees$p_crownHeight_mid <- NA

  new_trees$BA_mid<- NA
  new_trees$BAI_mid<- NA
  new_trees$weight_mid<- NA
  new_trees$height_mid<- NA
  new_trees$crownHeight_mid<- NA
  new_trees$volume_mid<- NA

  if (include_climate == FALSE){
    new_trees$t_avg <- NA
    new_trees$p_sum <- NA

  }

  new_trees$BAL_mid <- NA
  new_trees$stand_BA_mid <- NA
  new_trees$stand_n_mid <- NA

  protected_df <- df_before %>% dplyr::group_by(plotID) %>% dplyr::summarise(protected = median(protected))

  new_trees <- merge(new_trees, protected_df, by = "plotID", all.x = TRUE)

  new_trees <-  dplyr::select(new_trees, colnames(df_before))

  both_df <- rbind(df_before, new_trees)
  both_df$BAL <- NA
  both_df$stand_BA <- NA
  both_df$stand_n <- NA

  final_output_list <- list(

    predicted_ingrowth = both_df,
    eval_ingrowth = df_eval_ingrowth
  )

  # add the remaining elements
  for (i in 1:length(ing_codes)){

    final_output_list[[i+2]] <- get(paste0("mod_ing_", ing_codes[i]))
    names(final_output_list)[[i+2]] <- paste0("mod_ingrowth_",ing_codes[i])
  }

  return(final_output_list)

}
