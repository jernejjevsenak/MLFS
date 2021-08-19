#' volume_whole_tree_halfPeriod
#'
#' two-parameter volume functions for the MLFS
#' @keywords internal


volume_whole_tree_halfPeriod <- function(df) {

  # Define global variables
  BA_mid <- NULL


  initial_colnames <- colnames(df)

  df <- mutate(df,
               D = sqrt(4*BA_mid/pi) * 100
  )

  df$volume_mid <-

    ifelse(df[, "species"] == "PISY", #rdeci bor
           ( 0.35862280e+00*df[, "D"]*df[, "height"] - 0.57426745e-01*df[, "D"]*df[, "height"]^2 + 0.17327556e-02*df[, "D"]*df[, "height"]^3 + 0.68339275e-01*df[, "D"]^2*df[, "height"] - 0.10565319e-02*df[, "D"]^2*df[, "height"]^2 - 0.9992395e-04*df[, "D"]^3*df[, "height"] - 0.2057596e-06*df[, "D"]^3*df[, "height"]^3 + 0.17315704e-06*df[, "D"]^4*df[, "height"]^2 - 0.12513104e-10*df[, "D"]^5*df[, "height"]^3)/1000,

           ifelse(df[, "species"] == "LADE",  #macesen
                  ( 0.25691220e+00*df[, "D"]*df[, "height"] - 0.57880798e-02*df[, "D"]*df[, "height"]^2 + 0.26574631e-03*df[, "D"]*df[, "height"]^3 + 0.37477809e-01*df[, "D"]^2*df[, "height"] + 0.12434207e-03*df[, "D"]^2*df[, "height"]^2 - 0.31731646e-03*df[, "D"]^3*df[, "height"] - 0.144476078e-07*df[, "D"]^3*df[, "height"]^3 + 0.57839279e-07*df[, "D"]^4*df[, "height"]^2 - 0.33649650e-11*df[, "D"]^5*df[, "height"]^3)/1000,

                  ifelse(df[, "species"] == "PINI", #crni bor
                         ( 0.34301600e+00*df[, "D"]*df[, "height"] - 0.63656863e-01*df[, "D"]*df[, "height"]^2 + 0.20609616e-02*df[, "D"]*df[, "height"]^3 + 0.78368858e-01*df[, "D"]^2*df[, "height"] - 0.11691879e-02*df[, "D"]^2*df[, "height"]^2 - 0.51824098e-03*df[, "D"]^3*df[, "height"] - 0.68647134e-07*df[, "D"]^3*df[, "height"]^3 + 0.38480862e-06*df[, "D"]^4*df[, "height"]^2 - 0.85183903e-10*df[, "D"]^5*df[, "height"]^3)/1000,

                         ifelse(df[, "species"] == "ALGL" & df[, "D"] > 80, #jelsa
                                (pi * (df[, "D"]/20)^2 * (df[, "height"]*10)/2/1000 ), #jelsa valj
                                ifelse(df[, "species"] == "ALGL", #jelsa
                                       ( 0.20890280e-02*df[, "D"]*df[, "height"] + 0.89040978e-02*df[, "D"]*df[, "height"]^2 - 0.12841960e-03*df[, "D"]*df[, "height"]^3 + 0.28189206e-01*df[, "D"]^2*df[, "height"] - 0.36378839e-03*df[, "D"]^2*df[, "height"]^2 + 0.93421703e-03*df[, "D"]^3*df[, "height"] + 0.31139274e-06*df[, "D"]^3*df[, "height"]^3 - 0.10645016e-05*df[, "D"]^4*df[, "height"]^2 + 0.35466792e-09*df[, "D"]^5*df[, "height"]^3)/1000,

                                       ifelse(df[, "species"] == "FREX", #jesen (kot za bukev ker ni drugega)
                                              ( 0.17521850e+00*df[, "D"]*df[, "height"] - 0.20924022e-01*df[, "D"]*df[, "height"]^2 + 0.34280985e-03*df[, "D"]*df[, "height"]^3 + 0.58386328e-01*df[, "D"]^2*df[, "height"] - 0.27942794e-03*df[, "D"]^2*df[, "height"]^2 - 0.12736530e-03*df[, "D"]^3*df[, "height"] + 0.23287694e-07*df[, "D"]^3*df[, "height"]^3 + 0.88879276e-07*df[, "D"]^4*df[, "height"]^2 - 0.18411608e-10*df[, "D"]^5*df[, "height"]^3)/1000,

                                              ifelse(df[, "speciesGroup"] == "BEPE", #breza
                                                     ( 0.10108740e+00*df[, "D"]*df[, "height"] - 0.13010924e-01*df[, "D"]*df[, "height"]^2 + 0.59364890e-05*df[, "D"]*df[, "height"]^3 + 0.75920994e-01*df[, "D"]^2*df[, "height"] - 0.89132702e-03*df[, "D"]^2*df[, "height"]^2 - 0.11755791e-02*df[, "D"]^3*df[, "height"] - 0.22851593e-06*df[, "D"]^3*df[, "height"]^3 + 0.11332917e-05*df[, "D"]^4*df[, "height"]^2 - 0.28567640e-09*df[, "D"]^5*df[, "height"]^3)/1000,

                                                     ifelse(df[, "speciesGroup"] == 1,  #smreka
                                                            ( 0.2452791e+00*df[, "D"]*df[, "height"] - 0.1609334e-01*df[, "D"]*df[, "height"]^2 + 0.25960567e-03*df[, "D"]*df[, "height"]^3 + 0.55778174e-01*df[, "D"]^2*df[, "height"] - 0.24898613e-03*df[, "D"]^2*df[, "height"]^2 + 0.22516935e-03*df[, "D"]^3*df[, "height"] - 0.10924494e-07*df[, "D"]^3*df[, "height"]^3 + 0.36412393e-07*df[, "D"]^4*df[, "height"]^2 - 0.61708995e-12*df[, "D"]^5*df[, "height"]^3)/1000,

                                                            ifelse(df[, "speciesGroup"] == 2, #jelka
                                                                   ( 0.19819580e+00*df[, "D"]*df[, "height"] - 0.16444325e-01*df[, "D"]*df[, "height"]^2 + 0.25488555e-03*df[, "D"]*df[, "height"]^3 + 0.57219088e-01*df[, "D"]^2*df[, "height"] - 0.33008640e-03*df[, "D"]^2*df[, "height"]^2 + 0.18112376e-04*df[, "D"]^3*df[, "height"] - 0.77266170e-07*df[, "D"]^3*df[, "height"]^3 + 0.93429794e-08*df[, "D"]^4*df[, "height"]^2 - 0.50060347e-12*df[, "D"]^5*df[, "height"]^3)/1000,

                                                                   ifelse(df[, "speciesGroup"] == 4,  # bukev
                                                                          ( 0.17521850e+00*df[, "D"]*df[, "height"] - 0.20924022e-01*df[, "D"]*df[, "height"]^2 + 0.34280985e-03*df[, "D"]*df[, "height"]^3 + 0.58386328e-01*df[, "D"]^2*df[, "height"] - 0.27942794e-03*df[, "D"]^2*df[, "height"]^2 - 0.12736530e-03*df[, "D"]^3*df[, "height"] + 0.23287694e-07*df[, "D"]^3*df[, "height"]^3 + 0.88879276e-07*df[, "D"]^4*df[, "height"]^2 - 0.18411608e-10*df[, "D"]^5*df[, "height"]^3)/1000,

                                                                          ifelse(df[, "speciesGroup"] == 5, #hrasti
                                                                                 ( 0.17622750e+00*df[, "D"]*df[, "height"] - 0.20495698e-01*df[, "D"]*df[, "height"]^2 + 0.73239496e-03*df[, "D"]*df[, "height"]^3 + 0.53284931e-01*df[, "D"]^2*df[, "height"] - 0.92864308e-03*df[, "D"]^2*df[, "height"]^2 + 0.30671801e-04*df[, "D"]^3*df[, "height"] + 0.92752705e-07*df[, "D"]^3*df[, "height"]^3 - 0.80544290e-07*df[, "D"]^4*df[, "height"]^2 + 0.57597714e-11*df[, "D"]^5*df[, "height"]^3)/1000,



                                                                                 (pi * (df[, "D"]/20)^2 * (df[, "height"]*10)/2/1000 )))))))))))) #valj



  df <- select(df, all_of(initial_colnames))

  return(df)

}
