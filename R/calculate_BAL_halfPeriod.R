#' calculate_BAL_halfPeriod
#'
#' Function for the calculation of competition index BAL (Basal area in larger
#' trees)
#'
#' @keywords internal
#'

calculate_BAL_halfPeriod <- function(df){

  # Define global variables
  year <- NULL
  plotID <- NULL
  code <- NULL
  weight <- NULL
  weight_mid <- NULL
  BA <- NULL
  BA_mid <- NULL
  treeID <- NULL
  BA_ha <- NULL
  BA_ha_mid <- NULL
  count <- NULL
  BAL <- NULL
  BAL_mid <- NULL

  df$BAL_mid <- NA
  initial_colnames <- colnames(df)

  df$BAL_mid <- NULL

  # harvested trees get reduced weight
  temp <- mutate(df,
                 weight_mid = ifelse(code %in% c(1), weight_mid /2, weight_mid),
                 BA_ha_mid = BA_mid * weight_mid)

  temp <- select(temp, year, plotID, treeID, BA_ha_mid)

  temp <- temp %>% group_by(year, plotID) %>% mutate(count = row_number(plotID)) # %>% arrange(year, plotID, count)

  temp_sum <- reshape2::dcast(data = temp, formula = year + plotID ~ count, value.var = "BA_ha_mid")

  joined <- merge(temp, temp_sum, by = c("year", "plotID"))

  joined_BAL <- select(joined, -year, -plotID, -treeID, -count)

  joined_BAL[,-1][is.na(joined_BAL[,-1])] <- 0

  joined_BAL$BAL_mid <- rowSums(joined_BAL[-1] * (joined_BAL[,-1] >= joined_BAL[,1]), na.rm = TRUE)

  joined_BAL <- mutate(joined_BAL, BAL_mid = BAL_mid - BA_ha_mid)

  joined$BAL_mid <- joined_BAL$BAL_mid

  # final <- cbind(joined, joined_BAL[,"BAL"])

  final <- select(joined, year, plotID, treeID, BAL_mid)

  # summary(final)

  df1 <- merge(df, final, by = c("year", "plotID", "treeID"))

  df1 <- select(df1, all_of(initial_colnames))

  colnames(df1)

  return(df1)

}
