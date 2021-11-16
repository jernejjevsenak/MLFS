#' vol_form_factors
#'
#' two-parameter volume functions for the MLFS
#' @keywords internal

vol_form_factors <- function(df, form_factors = NULL, form_factors_level = "species",
                             uniform_form_factor = 0.42){

  initial_colnames <- colnames(df)

  # In case of form_factors_level = "species", make sure you have species in form factors data frame
  if (form_factors_level == "species"){

    if (!('species' %in% colnames(form_factors))){
      stop("column 'species' is missing in form_factors data frame")
    }

    if ('plotID' %in% colnames(form_factors)){
      stop("column 'plotID' should not be present in form_factors data frame")
    }
  }

  # In case of form_factors_level = "species_plot", make sure you have species and plotID in form factors data frame
  if (form_factors_level == "species_plot"){

    if (!('species' %in% colnames(form_factors))){
      stop("column 'species' is missing in form_factors data frame")
    }

    if (!('plotID' %in% colnames(form_factors))){
      stop("column 'plotID' is missing in form_factors data frame")
    }
  }

  # specify form factors
  if (is.null(form_factors)){

    df$form <- uniform_form_factor

  } else if (form_factors_level == "species"){

      if (!('species' %in% colnames(form_factors))){

        stop("form_factors should contain column 'species'")

      }

      df <- merge(df, form_factors, by = "species") # You can also add here third option, merge by plot & species

    } else if (form_factors_level == "plot"){

        if (!('plotID' %in% colnames(form_factors))){

          stop("form_factors should contain column 'plotID'")

        }


      df <- merge(df, form_factors, by = "plotID") # You can also add here third option, merge by plot & species

    } else if (form_factors_level == "species_plot"){

        if (sum(c('plotID','species') %in% colnames(form_factors)) < 1){

          stop("form_factors should contain columns 'plotID' and 'species'")

        }

      df <- merge(df, form_factors, by = c("plotID", "species")) # You can also add here third option, merge by plot & species

    } else {

      stop(paste0("form_factors_level should be 'species', 'plot' or 'species_plot' but instead it is ",  form_factors_level))

  }

  df$volume <- df$height * df$BA * df$form
  df$p_volume <- df$p_height * df$p_BA * df$form

  df <- dplyr::select(df, all_of(initial_colnames))

  return(df)
}
