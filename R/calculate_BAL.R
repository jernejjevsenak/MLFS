#' calculate_BAL
#'
#' Function for the calculation of competition index BAL (Basal area in larger
#' trees)
#'
#' @keywords internal
#'

calculate_BAL <- function(df){

  # Define global variables
  year <- NULL
  plotID <- NULL
  code <- NULL
  weight <- NULL
  BA <- NULL
  treeID <- NULL
  BA_ha <- NULL
  count <- NULL
  BAL <- NULL

  initial_colnames <- colnames(df)

  df$BAL <- NULL

  # harvested trees get reduced weight
  temp <- mutate(df, weight = ifelse(code %in% c(1), weight /2, weight),
                     BA_ha = BA * weight)

  temp <- select(temp, year, plotID, treeID, BA_ha)

  temp <- temp %>% group_by(year, plotID) %>% mutate(count = row_number(plotID)) # %>% arrange(year, plotID, count)

  temp_sum <- reshape2::dcast(data = temp, formula = year + plotID ~ count, value.var = "BA_ha")

  joined <- merge(temp, temp_sum, by = c("year", "plotID"))

  joined_BAL <- select(joined, -year, -plotID, -treeID, -count)

  joined_BAL[,-1][is.na(joined_BAL[,-1])] <- 0

  joined_BAL$BAL <- rowSums(joined_BAL[-1] * (joined_BAL[,-1] >= joined_BAL[,1]), na.rm = TRUE)

  joined_BAL <- mutate(joined_BAL, BAL = BAL - BA_ha)

  joined$BAL <- joined_BAL$BAL

  # final <- cbind(joined, joined_BAL[,"BAL"])

  final <- select(joined, year, plotID, treeID, BAL)

  # summary(final)

  df1 <- merge(df, final, by = c("year", "plotID", "treeID"))

  df1 <- select(df1, all_of(initial_colnames))

  colnames(df1)

  return(df1)

}
