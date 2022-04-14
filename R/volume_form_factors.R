#' volume_form_factors
#'
#' The calculation of individual tree volume using form factors, which can be
#' defined per species, per plot, or per species and per plot
#'
#' @param df data frame with tree heights and basal areas for individual trees
#' @param form_factors data frame with for factors for species, plot or both
#' @param form_factors_level character, the level of specified form factors. It
#' can be 'species', 'plot' or 'species_plot'
#' @param uniform_form_factor a uniform form factor to be applied to all trees.
#' If specified, it overwrites the argument 'form_factors'
#'
#' @return a data frame with calculated volume for all trees
#'
#' @examples
#' library(MLFS)
#' data(data_v3)
#' data(form_factors)
#'
#' data_v3 <- volume_form_factors(df = data_v3, form_factors = form_factors,
#'   form_factors_level = "species_plot")
#'
#' summary(data_v3)
#'

volume_form_factors <- function(df, form_factors = NULL, form_factors_level = "species",
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

      df <- merge(df, form_factors, by = "species")

    } else if (form_factors_level == "plot"){

        if (!('plotID' %in% colnames(form_factors))){

          stop("form_factors should contain column 'plotID'")

        }


      df <- merge(df, form_factors, by = "plotID")

    } else if (form_factors_level == "species_plot"){

        if (sum(c('plotID','species') %in% colnames(form_factors)) < 1){

          stop("form_factors should contain columns 'plotID' and 'species'")

        }

      df <- merge(df, form_factors, by = c("plotID", "species"))
    } else {

      stop(paste0("form_factors_level should be 'species', 'plot' or 'species_plot' but instead it is ",  form_factors_level))

  }

  df$volume <- df$height * df$BA * df$form
  df$p_volume <- df$p_height * df$p_BA * df$form

  df <- dplyr::select(df, all_of(initial_colnames))

  return(df)
}
