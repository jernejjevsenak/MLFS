#' vol_form_factors_halfPeriod
#'
#' two-parameter volume functions for the MLFS
#' @keywords internal

vol_form_factors_halfPeriod <- function(df, form_factors = NULL, form_factors_level = "species",
                             uniform_form_factor = 0.42){

  volume_mid <- NULL

  initial_colnames <- colnames(df)

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

  df$volume_mid <- df$height_mid * df$BA_mid * df$form

  df <- select(df, all_of(initial_colnames), volume_mid)

  return(df)
}
