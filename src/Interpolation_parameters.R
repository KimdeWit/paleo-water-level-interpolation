# Kim de Wit
# The code  presented in this script is part of the 
# “Living on Soft Soils: Subsidence and Society” project 
# (grant no. NWA.1160.18.259), WP1.3

###########################################################################
# Interpolation parameters
###########################################################################
# This script contains interpolation parameters fitted to 
# the HOLSEA-NL data set (v.3) and used in further interpolations

################################ Libraries ################################
library(gstat)

########################## Trend parameters ###############################


trend_outputs <- list(
  params_GW_withBgTect = c(a0 = 23.58352, a1 = 2.41678, a2 = -11.47674,
                           b = 2.95861,
                           c0 = 0.56955, c1 = -0.18611, c2 = -0.49841),
  
  params_GW_BgTectRem = c(a0 = 46.3201155, a1 = -15.6050644, a2 = -8.9359314,
                          b = 3.2571001,
                          c0 = 0.4653167, c1 = 0.3653018, c2 = -0.7266846),
  
  params_SL_withBgTect = c(a0 = 12.8543, a1 = 10.4087, a2 = -13.2684,
                           b = 2.6144,
                           c0 = 0.7000, c1 = -2.3303, c2 = 1.1914),
  
  params_SL_BgTectRem = c(a0 = 21.43464, a1 = 6.84338, a2 = -15.88182,
                          b =  2.77338,
                          c0 = 0.63654, c1 = -1.99606, c2 = 1.18340)
  
)


########################## Variograms #####################################

variograms <- list(

  varm_GW_withBgTect = vgm(psill = 0.94334, model = "Sph", 
                           nugget = 0.0599031,
                           range = 0.25, 
                           anis = c(90, 90, 90, 0.8, 0.5)),
                   
  varm_GW_BgTectRem = vgm(psill = 0.8617167, model = "Sph", 
                          nugget = 0.09131652,
                          range = 0.25, 
                          anis = c(90, 90, 90, 0.8, 0.5)),
  
  varm_SL_withBgTect = vgm(psill = 0.8709085, model = "Sph",
                           nugget =  0.09499782,
                           range = 0.20,
                           anis = c(90, 90, 90, 0.8, 0.5)),
                   
  varm_SL_BgTectRem = vgm(psill = 0.7124101, model = "Sph", 
                          nugget = 0.1017607,
                          range = 0.2, 
                          anis = c(90, 90, 90, 0.8, 0.5))
)



# with limits variogram results:
# psill
# SL-BgTectRem = 0.7754955

# SL-withBgTect = 1.14799

# GW-BgTectRem = 0.8607831

# GW-withBgTect = 0.9512526
