#' V_general_halfPeriod
#'
#' three-parameter volume functions for the MLFS
#' @keywords internal

V_general_halfPeriod = function(df, data_volF_param = data_volF_param){

  species <- NULL
  BA_mid <- NULL
  p_BA <- NULL
  height_mid <- NULL
  p_height <- NULL
  H <- NULL
  p_D <- NULL
  p_H <- NULL
  plotID <- NULL
  treeID <- NULL
  equation <- NULL
  equation_mid <- NULL
  rowwise <- NULL
  DBH_min <- NULL
  DBH_max <- NULL
  vol_factor <- NULL
  BA <- NULL
  BA_min <- NULL
  BA_max <- NULL
  valid <- NULL

  initial_colnames <- colnames(df)

  # I create new variable 'species temp', which is used to merge with volume functions
  df <- dplyr::mutate(df, species_temp = ifelse(species %in% unique(data_volF_param$species), species, "REST"))
  data_volF_param <- dplyr::rename(data_volF_param, 'species_temp' = 'species')

  # Check that we have all equations or, that we have the REST equation available
  if (sum(!(unique(df$species) %in% unique(data_volF_param$species_temp))) > 1){

    if (!("REST" %in% unique(data_volF_param$species_temp))){

      stop(paste0("The equations and parameters for volume functions are missing for the following species:",
                  paste0(unique(df$species)[!(unique(df$species) %in% unique(data_volF_param$species_temp))], collapse=", "),
                  ". Alternatively, provide the equation for the 'REST' category"))
    }
  }

  # If volume functions are provided per plot and species, we do merge using both variables, otherwise only species (species_temp)
  if ("plotID" %in% colnames(data_volF_param)){

    df <- merge(df, data_volF_param, by = c("plotID", "species_temp"), all.x = TRUE)

  } else {

    df <- merge(df, data_volF_param, by = "species_temp", all.x = TRUE)

  }

  df <- df %>%
    rowwise() %>%

      mutate(

      BA_min = ((DBH_min/2)^2 * pi)/10000,
      BA_max = ((DBH_max/2)^2 * pi)/10000,

      D_mid = sqrt(4*BA_mid/pi) * 100, #cm
      H_mid = height_mid, # m

      equation_mid = gsub("H", "H_mid", gsub("D", "D_mid", equation)),
      volume_mid = eval(parse(text=equation_mid))/vol_factor,

      valid = ifelse(is.na(DBH_min), TRUE,
                     ifelse(BA >= BA_min &  BA <= BA_max, TRUE, FALSE))
      ) %>% filter(valid == TRUE)

  df <- dplyr::select(df, all_of(initial_colnames)) %>% arrange(plotID, treeID)

  df <- data.frame(df)

  return(df)

}
