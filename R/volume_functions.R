#' V_general
#'
#' three-parameter volume functions for the MLFS
#' @keywords internal

V_general = function(df, data_volF_param = data_volF_param){

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

  initial_colnames <- colnames(df)

  #####################
  # 1 Volume for PCAB #
  #####################
  df_PCAB <- filter(df, species == "PCAB")
  param_PCAB <- filter(data_volF_param, species == "PCAB")

  a = param_PCAB$a
  b = param_PCAB$b
  c = param_PCAB$c
  d = param_PCAB$d
  e = param_PCAB$e
  f = param_PCAB$f
  g = param_PCAB$g

  df_PCAB <- mutate(df_PCAB,
               D = sqrt(4*BA/pi) * 100/10, #cm -> dm
               p_D = sqrt(4*p_BA/pi) * 100/10, #cm -> dm
               H = height * 10, # m -> dm
               p_H = p_height * 10,
               volume = (pi/4)*(a*D^2*H+b*D^2*H*log(D)^2+c*D^2+d*D*H+e*H+f*D) / 1000,   # dm3 -> m3
               volume =  (pi/4)*(a*D^2*H+b*D^2*H*log(D)^2+c*D^2+d*D*H+e*H+f*D+g)/1000,   # dm3 -> m3
               p_volume = (pi/4)*(a*p_D^2*p_H+b*p_D^2*p_H*log(p_D)^2+c*p_D^2+d*p_D*p_H+e*p_H+f*D) / 1000)   # dm3 -> m3

  #####################
  # 2 Volume for ABAL #
  #####################
  df_ABAL <- filter(df, species == "ABAL")
  param_ABAL <- filter(data_volF_param, species == "ABAL")

  a = param_ABAL$a
  b = param_ABAL$b
  c = param_ABAL$c
  d = param_ABAL$d
  e = param_ABAL$e
  f = param_ABAL$f
  g = param_ABAL$g

  df_ABAL <- mutate(df_ABAL,
                    D = sqrt(4*BA/pi) * 100/10, # cm -> dm
                    p_D = sqrt(4*p_BA/pi) * 100/10, # cm -> dm
                    H = height * 10, # m -> dm
                    p_H = p_height * 10,
                    volume = (pi/4)*(a*D^2*H+b*D^2*H*log(D)^2+c*D^2+d*D*H+e*H+f*D+g)/1000,  # dm3 -> m3
                    p_volume = (pi/4)*(a*p_D^2*p_H+b*p_D^2*p_H*log(p_D)^2+c*p_D^2+d*p_D*p_H+e*p_H+f*p_D+g)/1000)   # dm3 -> m3


  #####################
  # 3 Volume for FASY #
  #####################
  df_FASY <- filter(df, species == "FASY")
  param_FASY <- filter(data_volF_param, species == "FASY")

  a = param_FASY$a
  b = param_FASY$b
  c = param_FASY$c
  d = param_FASY$d
  e = param_FASY$e
  f = param_FASY$f
  g = param_FASY$g

  df_FASY <- mutate(df_FASY,
                    D = sqrt(4*BA/pi) * 100/10, # cm -> dm
                    p_D = sqrt(4*p_BA/pi) * 100/10, # cm -> dm
                    H = height * 10, # m -> dm
                    p_H = p_height * 10,
                    volume = (pi/4)*(a*D^2*H+b*D^2*H*log(D)^2+c*D^2+d*D*H+e*H+f*D+g)/1000,   # dm3 -> m3
                    p_volume = (pi/4)*(a*p_D^2*p_H+b*p_D^2*p_H*log(p_D)^2+c*p_D^2+d*p_D*p_H+e*p_H+f*p_D+g)/1000)   # dm3 -> m3


  #####################
  # 4 Volume for PISY #
  #####################
  df_PISY <- filter(df, species == "PISY")
  param_PISY <- filter(data_volF_param, species == "PISY")

  a = param_PISY$a
  b = param_PISY$b
  c = param_PISY$c
  d = param_PISY$d
  e = param_PISY$e
  f = param_PISY$f
  g = param_PISY$g

  df_PISY <- mutate(df_PISY,
                    D = sqrt(4*BA/pi) * 100/10, # cm -> dm
                    p_D = sqrt(4*p_BA/pi) * 100/10, # cm -> dm
                    H = height * 10, # m -> dm
                    p_H = p_height * 10,
                    volume = (pi/4)*(a*D^2*H+b*D^2*H*log(D)^2+c*D^2+d*D*H)/1000,   # dm3 -> m3
                    p_volume = (pi/4)*(a*p_D^2*p_H+b*p_D^2*p_H*log(p_D)^2+c*p_D^2+d*p_D*p_H)/1000)   # dm3 -> m3


  #####################
  # 5 Volume for QUSP #
  #####################
  df_QUSP <- filter(df, species %in% c("QUSP", "QUCE", "QUPE", "QUPU", "QURO", "QURU", "QUSP"))
  param_QUSP <- filter(data_volF_param, species == "QUSP")

  a = param_QUSP$a
  b = param_QUSP$b
  c = param_QUSP$c
  d = param_QUSP$d
  e = param_QUSP$e
  f = param_QUSP$f
  g = param_QUSP$g

  df_QUSP <- mutate(df_QUSP,
                    D = sqrt(4*BA/pi) * 100/10, #cm -> dm
                    p_D = sqrt(4*p_BA/pi) * 100/10, #cm -> dm
                    H = height * 10, # m -> dm
                    p_H = p_height * 10,
                    volume = (pi/4)*(a*D^2*H+b*D^2+c*D*H+d*H+e*D+f) / 1000,   # dm3 -> m3
                    p_volume = (pi/4)*(a*p_D^2*p_H+b*p_D^2+c*p_D*p_H+d*p_H+e*p_D+f) / 1000)   # dm3 -> m3

  #####################
  # 4 Volume for Rest #
  #####################

  df_REST <- filter(df, !(species %in% c("PCAB", "ABAL", "FASY", "PISY",
                                         "QUSP", "QUCE", "QUPE", "QUPU", "QURO", "QURU", "QUSP")))

  param_REST <- filter(data_volF_param, species == "REST") # is this the best option?

  a = param_REST$a
  b = param_REST$b
  c = param_REST$c
  d = param_REST$d
  e = param_REST$e
  f = param_REST$f
  g = param_REST$g

  df_REST <- mutate(df_REST,
                    D = sqrt(4*BA/pi) * 100/10, #cm -> dm
                    p_D = sqrt(4*p_BA/pi) * 100/10, #cm -> dm
                    H = height * 10, # m -> dm
                    p_H = p_height * 10,
                    volume = (pi/4)*(a*D^2*H+b*D^2*H*log(D)^2+c*D^2+d*D*H+e*H+f*D) / 1000,   # dm3 -> m3
                    p_volume = (pi/4)*(a*p_D^2*p_H+b*p_D^2*p_H*log(p_D)^2+c*p_D^2+d*p_D*p_H+e*p_H+f*D) / 1000)   # dm3 -> m3

  ########
  # Join #
  ########

  df <- rbind(df_PCAB,
              df_ABAL,
              df_FASY,
              df_PISY,
              df_QUSP,
              df_REST)

  df <- select(df, all_of(initial_colnames)) %>% arrange(plotID, treeID)

  return(df)
}
