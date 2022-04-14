#' volume_merchantable_halfPeriod
#'
#' two-parameter volume functions for the MLFS
#'
#' @return a data frame with calculated volume for all trees in the middle of
#' a simulation step
#'
#' @keywords internal

volume_merchantable_halfPeriod <- function(df) {

  # Define global variables
  BA_mid <- NULL

  initial_colnames <- colnames(df)

  df <- mutate(df, D = sqrt(4*BA_mid/pi) * 100)

  df$volume_mid <- ifelse(df[, "species"] == "PISY", #rdeci bor
           (-0.25718650e+00*df[, "D"]*df[, "height"] + 0.83744137e-02*df[, "D"]*df[, "height"]^2 + 0.48990826e-04*df[, "D"]*df[, "height"]^3 + 0.53526983e-01*df[, "D"]^2*df[, "height"] - 0.71004612e-03*df[, "D"]^2*df[, "height"]^2 - 0.64284461e-04*df[, "D"]^3*df[, "height"] + 0.54395432e-07*df[, "D"]^3*df[, "height"]^3 + 0.78186279e-07*df[, "D"]^4*df[, "height"]^2 - 0.21088095e-10*df[, "D"]^5*df[, "height"]^3)/1000,

           ifelse(df[, "species"] == "PINI", #crni bor
                  (-0.14904420e+00*df[, "D"]*df[, "height"] - 0.11859180e-02*df[, "D"]*df[, "height"]^2 + 0.34528025e-03*df[, "D"]*df[, "height"]^3 + 0.57969716e-01*df[, "D"]^2*df[, "height"] - 0.71422313e-03*df[, "D"]^2*df[, "height"]^2 - 0.29259742e-03*df[, "D"]^3*df[, "height"] + 0.52762381e-06*df[, "D"]^3*df[, "height"]^3 + 0.16490879e-06*df[, "D"]^4*df[, "height"]^2 - 0.29892899e-10*df[, "D"]^5*df[, "height"]^3)/1000,

                  ifelse(df[, "species"] == "LADE",  #macesen
                         (-0.11595690e+00*df[, "D"]*df[, "height"] + 0.17085340e-01*df[, "D"]*df[, "height"]^2 - 0.51256211e-04*df[, "D"]*df[, "height"]^3 + 0.32234182e-01*df[, "D"]^2*df[, "height"] + 0.13366208e-03*df[, "D"]^2*df[, "height"]^2 - 0.31517344e-03*df[, "D"]^3*df[, "height"] - 0.19920908e-07*df[, "D"]^3*df[, "height"]^3 + 0.57817373e-07*df[, "D"]^4*df[, "height"]^2 - 0.45907137e-11*df[, "D"]^5*df[, "height"]^3)/1000,

                         ifelse(df[, "species"] == "ALGL" & df[, "D"] > 80, #jelsa
                                (pi * (df[, "D"]/20)^2 * (df[, "height"]*10)/2/1000 ), #jelsa valj
                                ifelse(df[, "species"] == "ALGL", #jelsa
                                       (-0.55012480e+00*df[, "D"]*df[, "height"] + 0.29872117e-01*df[, "D"]*df[, "height"]^2 - 0.21623699e-03*df[, "D"]*df[, "height"]^3 + 0.76450451e-01*df[, "D"]^2*df[, "height"] - 0.16553548e-02*df[, "D"]^2*df[, "height"]^2 - 0.44600282e-03*df[, "D"]^3*df[, "height"] + 0.48874635e-06*df[, "D"]^3*df[, "height"]^3 + 0.29245893e-06*df[, "D"]^4*df[, "height"]^2 - 0.13995810e-09*df[, "D"]^5*df[, "height"]^3)/1000,

                                       ifelse(df[, "speciesGroup"] == "BEPE", #breza
                                              (-0.48723660e+00*df[, "D"]*df[, "height"] + 0.24600888e-01*df[, "D"]*df[, "height"]^2 + 0.53121062e-05*df[, "D"]*df[, "height"]^3 + 0.81455648e-01*df[, "D"]^2*df[, "height"] - 0.18627497e-02*df[, "D"]^2*df[, "height"]^2 - 0.10328362e-02*df[, "D"]^3*df[, "height"] + 0.40198766e-06*df[, "D"]^3*df[, "height"]^3 + 0.98021908e-06*df[, "D"]^4*df[, "height"]^2 - 0.35491548e-09*df[, "D"]^5*df[, "height"]^3)/1000,

                                              ifelse(df[, "species"] == "FREX", #jesen
                                                     (-0.53354550e-01*df[, "D"]*df[, "height"] - 0.23240470e-02*df[, "D"]*df[, "height"]^2 + 0.73296776e-04*df[, "D"]*df[, "height"]^3 + 0.42556037e-01*df[, "D"]^2*df[, "height"] + 0.76818087e-04*df[, "D"]^2*df[, "height"]^2 - 0.16196035e-03*df[, "D"]^3*df[, "height"] - 0.96165262e-07*df[, "D"]^3*df[, "height"]^3 + 0.30836228e-07*df[, "D"]^4*df[, "height"]^2 + 0.20132201e-10*df[, "D"]^5*df[, "height"]^3)/1000,

                                                     ifelse(df[, "speciesGroup"] == 1,  #smreka
    (-0.2395044e+00*df[, "D"]*df[, "height"] + 0.68976337e-02*df[, "D"]*df[, "height"]^2 - 0.37308357e-05*df[, "D"]*df[, "height"]^3 + 0.59380420e-01*df[, "D"]^2*df[, "height"] - 0.30222736e-03*df[, "D"]^2*df[, "height"]^2 - 0.43109400e-03*df[, "D"]^3*df[, "height"] + 0.64316383e-08*df[, "D"]^3*df[, "height"]^3 + 0.90518929e-07*df[, "D"]^4*df[, "height"]^2 - 0.90550376e-11*df[, "D"]^5*df[, "height"]^3)/1000,

                ifelse(df[, "speciesGroup"] == 2, #jelka
    (-0.22822880e+00*df[, "D"]*df[, "height"] + 0.19664395e-01*df[, "D"]*df[, "height"]^2 - 0.28041458e-03*df[, "D"]*df[, "height"]^3 + 0.44570362e-01*df[, "D"]^2*df[, "height"] - 0.14288404e-03*df[, "D"]^2*df[, "height"]^2 - 0.11572040e-03*df[, "D"]^3*df[, "height"] + 0.30996593e-07*df[, "D"]^3*df[, "height"]^3 + 0.37584951e-08*df[, "D"]^4*df[, "height"]^2 - 0.18066064e-11*df[, "D"]^5*df[, "height"]^3)/1000,

                ifelse(df[, "speciesGroup"] == 4,  # bukev
    (-0.19403830e+00*df[, "D"]*df[, "height"] + 0.58124245e-02*df[, "D"]*df[, "height"]^2 - 0.30116138e-05*df[, "D"]*df[, "height"]^3 + 0.46951914e-01*df[, "D"]^2*df[, "height"] - 0.22621894e-03*df[, "D"]^2*df[, "height"]^2 - 0.10182239e-03*df[, "D"]^3*df[, "height"] + 0.43435806e-07*df[, "D"]^3*df[, "height"]^3 + 0.84446836e-07*df[, "D"]^4*df[, "height"]^2 - 0.17618793e-10*df[, "D"]^5*df[, "height"]^3)/1000,

                ifelse(df[, "speciesGroup"] == 5, #hrasti
    (-0.26834540e+00*df[, "D"]*df[, "height"] + 0.49444128e-02*df[, "D"]*df[, "height"]^2 + 0.10771745e-03*df[, "D"]*df[, "height"]^3 + 0.57851059e-01*df[, "D"]^2*df[, "height"] - 0.57706566e-03*df[, "D"]^2*df[, "height"]^2 - 0.84271339e-04*df[, "D"]^3*df[, "height"] + 0.42123926e-07*df[, "D"]^3*df[, "height"]^3 + 0.59138330e-07*df[, "D"]^4*df[, "height"]^2 - 0.98238961e-11*df[, "D"]^5*df[, "height"]^3)/1000,


    (pi * (df[, "D"]/20)^2 * (df[, "height"]*10)/2/1000 )))))))))))) # valj


  df <- dplyr::select(df, all_of(initial_colnames))

  return(df)
}


