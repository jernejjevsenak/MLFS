#' volume_functions
#'
#' The calculation of individual tree volume using the n-parameter volume
#' functions for the MLFS
#'
#' @param df data frame with tree heights and basal areas for individual trees
#' @param df_volumeF_parameters data frame with equations and parameters for
#' n-parametric volume functions
#'
#' @return a data frame with calculated volume for all trees
#'
#' @examples
#' library(MLFS)
#' data(data_v3)
#' data(df_volume_parameters)
#'
#' data_v3 <- volume_functions(df = data_v3,
#'   df_volumeF_parameters = df_volume_parameters)
#'

volume_functions = function(df, df_volumeF_parameters = NULL){

  species <- NULL
  BA <- NULL
  p_BA <- NULL
  height <- NULL
  p_height <- NULL
  H <- NULL
  p_D <- NULL
  p_H <- NULL
  plotID <- NULL
  treeID <- NULL
  equation <- NULL
  p_equation <- NULL
  rowwise <- NULL
  DBH_min <- NULL
  DBH_max <- NULL
  vol_factor <- NULL
  BA_min <- NULL
  BA_max <- NULL
  valid <- NULL
  p_volume <- NULL
  p_valid <- NULL
  d_factor <- NULL
  h_factor <- NULL

  initial_colnames <- colnames(df)

  # I create new variable 'species temp', which is used to merge with volume functions
  df <- dplyr::mutate(df, species_temp = ifelse(species %in% unique(df_volumeF_parameters$species), species, "REST"))
  df_volumeF_parameters <- dplyr::rename(df_volumeF_parameters, 'species_temp' = 'species')

  # Check that we have all equations or, that we have the REST equation available
  if (sum(!(unique(df$species) %in% unique(df_volumeF_parameters$species_temp))) > 1){

    if (!("REST" %in% unique(df_volumeF_parameters$species_temp))){

      stop(paste0("The equations and parameters for volume functions are missing for the following species:",
                  paste0(unique(df$species)[!(unique(df$species) %in% unique(df_volumeF_parameters$species_temp))], collapse=", "),
                  ". Alternatively, provide the equation for the 'REST' category"))
    }
  }

  # If volume functions are provided per plot and species, we do merge using both variables, otherwise only species (species_temp)
  if ("plotID" %in% colnames(df_volumeF_parameters)){

    df <- merge(df, df_volumeF_parameters, by = c("plotID", "species_temp"), all.x = TRUE)

  } else {

    df <- merge(df, df_volumeF_parameters, by = "species_temp", all.x = TRUE)

  }

  df_current <- df %>%

         rowwise() %>% mutate(

         BA_min = ((DBH_min/2)^2 * pi)/10000,
         BA_max = ((DBH_max/2)^2 * pi)/10000,

         D = sqrt(4*BA/pi) * 100/d_factor,
         H = height*h_factor,
         volume = eval(parse(text=equation))/vol_factor,

         valid = ifelse(is.na(DBH_min), TRUE,
                        ifelse(BA >= BA_min &  BA <= BA_max, TRUE, FALSE))

         ) %>% filter(valid == TRUE) %>% select(-p_BA, -p_height, -p_volume)

df_previous <- df %>%
    rowwise() %>% mutate(

      BA_min = ((DBH_min/2)^2 * pi)/10000,
      BA_max = ((DBH_max/2)^2 * pi)/10000,

      p_D = sqrt(4*p_BA/pi) * 100/d_factor,
      p_H = p_height * h_factor,
      p_equation = gsub("H", "p_H", gsub("D", "p_D", equation)),
      p_volume = eval(parse(text=p_equation))/vol_factor,
      p_valid = ifelse(is.na(DBH_min), TRUE,
                       ifelse(p_BA >= BA_min &  p_BA <= BA_max, TRUE, FALSE))

    ) %>% filter(p_valid == TRUE) %>% select(plotID, treeID, p_BA, p_height, p_volume)

  df <- merge(df_current, df_previous, by = c("plotID", "treeID"))

  df <- dplyr::select(df, all_of(initial_colnames)) %>% arrange(plotID, treeID)

  df <- data.frame(df)

  return(df)
}
