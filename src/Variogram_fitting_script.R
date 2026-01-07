# Kim de Wit
# The code  presented in this script is part of the 
# “Living on Soft Soils: Subsidence and Society” project 
# (grant no. NWA.1160.18.259), WP1.3

###########################################################################
# Kriging variogram plotting script
###########################################################################
# This script is designed to determine kriging variograms suitable for a 3D 
# groundwater level or sea-level interpolation for the Netherlands using the 
# HOLSEA-NL data set as input.

# The input data file needs to be a table with preferably x, y and z coordinates 
# with RD new coordinate projection. 

###########################################################################
########## Clear out environment before running the script ################
rm(list = ls())
gc()

################# Set working directory and library path ##################
# Set working directory and library path ----------------------------------
# Store output in organised folders and files using current month+year (folder)
# and the current day+month+year (output-files)
folder <- format(Sys.time(), "%B_%Y")

if (file.exists(file.path("./results/output", folder, "Interpolation_output", 
                          paste0(Sys.Date())))) {
  cat("The folder already exists")
} else if (file.exists(file.path("./results/output/", folder))) {
  dir.create(file.path("./results/output/", folder, "Interpolation_output", 
                       paste0(Sys.Date())), recursive = TRUE)
} else {
  dir.create(file.path("./results/output/", folder, "Interpolation_output", 
                       paste0(Sys.Date())), recursive = TRUE)
}

output_loc <- paste("./results/output", folder, "Interpolation_output",
                    paste0(Sys.Date()), sep = "/")


################################ Libraries ################################
library(gstat)
library(snakecase)
library(readxl)
library(sp)
library(stars)
library(sf)
library(dplyr)
library(tidyverse)

library(tictoc) # only for time check

source("./src/general_functions.R")
source("./src/interpolation_parameters.R")


tic()

############################ Prepare input data ###########################

# Set coordinate projection
RDnew <- st_crs("EPSG:28992")

# Observation data - prepare for trend fitting
# Loading HOLSEA excel sheet directly, removing top headers and N/A and nd values
GW_indexpoints_data <- read_excel("data/raw-HOLSEA/HOLSEA-NL-DeWit-Cohen-2024-full-dataset-v3.xlsx",
                                  sheet = "Long-form", na = c("nd", "n/a"), skip = 2)

GW_indexpoints_data <- adjust_headers(GW_indexpoints_data)

# Specify input data, Region and type of trend equation ----------------------

# Define which geological data types should be included in the analysis
# And if the sea-level or gw-level elevation should be used
# obs_type: type of indicators used from observations
#   c(0) = SLIPs only, c(0,1) = SLIPs & tidal ULD, 
#   c(0,1,2,3) = incl. river gradient & all ULDs
# indicator_type: determines if the groundwater ("GW") or sea-level ("SL")
#   elevation is used
# indicator: Label used for file name ("GW-level" or "sea-level")
# coast: if TRUE, the interpolation only uses data from the coastal plain (e.g.
#   where HDEM < 1m)

# Un-comment the preferred selection
# # GW
# obs_type = c(0,1,2,3)
# indicator_type = "GW"
# indicator = "GW-level"
# coast = FALSE # Standard FALSE for GW-level data

# Sea-level
obs_type = c(0,1)
indicator_type = "SL" 
indicator = "sea-level" 
coast = TRUE # Standard TRUE for Sea-level data

# correction: determines if a correction for tectonic subsidence is used
#   withBgTect: when excluding tectonic correction
#   BgTectRem: to correct data for the background tectonic signal
# correction = "withBgTect" 
correction = "BgTectRem" 

# Determine if coefficients should be floored at a minimal value
#   with_limits: Coefficient a imited to minimal value of 3 and coefficient c 
#                limited to minimal value of 0.
#   no_limits:   using the original fitted coefficient values
# coef_lim = "no_limits" 
coef_lim = "with_limits"

# Select regions to include (or leave out) in trend fit: 
# For the HOLSEA-NL dataset the same regions as in De Wit & Cohen (2024) 
# are used:
# "Noord-Holland", "Rhine-Meuse delta inland", "Flevoland", "Zuid-Holland", 
# "Zeeland", "Waddenzee West", "Waddenzee East"
# Standard all regions are used
regions = c("Noord-Holland", "Rhine-Meuse delta inland", 
            "Flevoland", "Zuid-Holland", "Zeeland",
            "Waddenzee West", "Waddenzee East")


############################### Functions #################################
# Define trend equation 

# The trend fit with a & c linear dependent on x and y coordinates is the
# standard option used for the interpolation in the Netherlands.
# If an other trend function is desired, this section and the parameters in
# interpolation_parameters.R would have to be adjusted accordingly.

# a & c linear
# (b is constant, a and c are linearly dependent on x and y)
fitmodel_name = "Linearfit_ac"

calc_a <- function(params, x_n, y_n){
  params[1] + params[2]*x_n + params[3]*y_n
}

calc_c <- function(params, x_n, y_n){
  params[5] + params[6]*x_n + params[7]*y_n
}


# Interpolation functions
# Define sigmoid function
sigmoid_function = function(a, b, c, q, p, x, y) {
  (1 - c) * (1 - exp(-a * q * p^b)) + (c * p)
}

# Function to check for duplicate (x, y) coordinates and modify them
# Kriging does not allow for duplicate coordinates, therefor they need to be 
# adjusted before kriging (adding 0.1 to x-coordinate for duplicates)
remove_duplicates <- function(df, x_col = "x", y_col = "y") {
  for (i in 1:nrow(df)) {
    while (sum(df[[x_col]] == df[[x_col]][i] & df[[y_col]] == df[[y_col]][i]) > 1) {
      df[i, x_col] <- df[i, x_col] + 0.1  
    }
  }
  return(df)
}

############################## Define Constants ###########################
Dmean = 11.5  # mean thickness between LDEM and HDEM
Dmax = 27.1   # max thickness between LDEM and HDEM
A0 = -10800   # [cal years Before Present (BP)] Start time of studied period 
A1 = -1000    # [cal years BP] end time of studied period 
tmin = -A0
tmax = -A1
tstep = 200   # [years]

t = seq(from = A0, to = A1, by = tstep) # Time steps
p = (A0-t)/(A0-A1)                      # Normalized time steps

# Grid range - based on extent DEM - in RD new coordinates
xmin = -6000   # [m]
xmax = 280000  # [m]
ymin = 340000  # [m]
ymax = 645000  # [m]
xstep = 10000   # [m]
ystep = 10000   # [m]

# Grid - based on extent DEMs
grid2D = expand.grid(x = seq (from = xmin, to = xmax, by = xstep),
                     y = seq(from = ymin, to = ymax, by = ystep))
grid3D = expand.grid(x = seq (from = xmin, to = xmax, by = xstep),
                     y = seq(from = ymin, to = ymax, by = ystep),
                     age = seq(from = -A0, to = -A1, by = -tstep))

grid2D_points <- st_as_sf(grid2D, coords = c("x","y"), crs = RDnew)

# Interpolation trend and kriging parameters
params_key <- paste("params", indicator_type, correction, sep = "_")
params <- trend_outputs[[params_key]]

# Define bounding surfaces and Pleistocene surface
# Each DEM has a uncorrected (withBgTect) and corrected (BgTectRem) version,
# depending on if the background tectonic subsidence is removed or not.
PDEM_withBgTect = "data/processed/DEMs/PDEM_extended.tif"
PDEM_BgTectRem  = "data/processed/DEMs/PDEM_extended_BgTectCorrected.tif"
HDEM_withBgTect = "data/processed/DEMs/HDEM_extended_NL_v2.tif"
LDEM_withBgTect = "data/processed/DEMs/LDEM_extended_NL_v2.tif"
HDEM_BgTectRem  = "data/processed/DEMs/HDEM_extended_BgTectCorrected.tif"
LDEM_BgTectRem  = "data/processed/DEMs/LDEM_extended_BgTectCorrected.tif"
VLM_rate = "data/processed/DEMs/Mean_VLM_rate_extended.tif"


############################ Process input data ###########################
# Select observation data from HOLSEA data set based on input selection
obs <- 
  if (indicator_type == "SL") {
    
    GW_indexpoints_data %>%
      filter(reject %in% 0) %>%
      filter(type %in% obs_type) %>%
      select(ID = labidnr, 
             Name = sample_name,
             x = rd_x, 
             y = rd_y,
             Z = corrected_rsl_m_if_any,
             Z_mean = corrected_rsl_m_if_any,
             Z_2σ_min = corrected_rsl_2_σ_uncertainty_m_min,
             Z_2σ_plus = corrected_rsl_2_σ_uncertainty_m_plus,
             Z_BgTect = corrected_rsl_with_bg_tect_m_if_any,
             Z_BgTect_mean = corrected_rsl_with_bg_tect_m_if_any,
             Z_BgTect_2σ_min = withBG_tect_corrected_rsl_2_σ_uncertainty_m_min,
             Z_BgTect_2σ_plus = withBG_tect_corrected_rsl_2_σ_uncertainty_m_plus,
             tectonic_cor = tectonic_correction_m,
             tectonic_err = tectonic_correction_uncertainty_m,
             type, 
             age_cal_a_bp,
             age_2_σ_uncertainty_cal_a_plus,
             age_2_σ_uncertainty_cal_a_min,
             ox_cal_modelled_age_cal_a_bp,
             ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus,
             ox_cal_modelled_age_2_σ_uncertainty_cal_a_min,
             region_name
      ) %>%
      mutate(
        t = if_else(!is.na(ox_cal_modelled_age_cal_a_bp) , 
                    ox_cal_modelled_age_cal_a_bp, age_cal_a_bp)
      ) %>%
      mutate(
        Age_2σ_plus = if_else(!is.na(ox_cal_modelled_age_2_σ_uncertainty_cal_a_min) , 
                              ox_cal_modelled_age_2_σ_uncertainty_cal_a_min, 
                              age_2_σ_uncertainty_cal_a_min)
      ) %>%
      mutate(
        Age_2σ_min = if_else(!is.na(ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus) , 
                             ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus, 
                             age_2_σ_uncertainty_cal_a_plus)
      ) %>%
      mutate(
        Z = if_else(type > 0, Z - Z_2σ_min, Z)
      ) %>%
      mutate(
        Z_BgTect = if_else(type > 0, Z_BgTect - Z_BgTect_2σ_min, Z_BgTect)
      )
  } else {
    GW_indexpoints_data %>%
      filter(reject %in% 0) %>%
      filter(type %in% obs_type) %>%
      mutate(
        human_induced_subsidence_uncertainty_m = 
          if_else(!is.na(human_induced_subsidence_uncertainty_m) , 
                  human_induced_subsidence_uncertainty_m, 0),
        compaction_correction_uncertainty_m = 
          if_else(!is.na(compaction_correction_uncertainty_m) , 
                  compaction_correction_uncertainty_m, 0)
      ) %>%
      mutate(Z_BgTect = corrected_z_gw_m - tectonic_correction_m,
             Z_BgTect_2σ_min = sqrt(sample_elevation_uncertainty_m_plus^2 + 
                                      indicative_range_uncertainty_m^2 +
                                      compaction_correction_uncertainty_m^2
                                    + human_induced_subsidence_uncertainty_m^2),
             Z_BgTect_2σ_plus = sqrt(sample_elevation_uncertainty_m_min^2 + 
                                       indicative_range_uncertainty_m^2 +
                                       compaction_correction_uncertainty_m^2
                                     + human_induced_subsidence_uncertainty_m^2)
      ) %>%
      select(ID = labidnr, 
             Name = sample_name,
             x = rd_x, 
             y = rd_y,
             Z = corrected_z_gw_m,
             Z_mean = corrected_z_gw_m,
             Z_2σ_min = corrected_gwl_2_σ_uncertainty_m_min,
             Z_2σ_plus = corrected_gwl_2_σ_uncertainty_m_plus,
             Z_BgTect = Z_BgTect,
             Z_BgTect_mean = Z_BgTect,
             Z_BgTect_2σ_min,
             Z_BgTect_2σ_plus,
             tectonic_cor = tectonic_correction_m,
             tectonic_err = tectonic_correction_uncertainty_m,
             type, 
             age_cal_a_bp,
             age_2_σ_uncertainty_cal_a_plus,
             age_2_σ_uncertainty_cal_a_min,
             ox_cal_modelled_age_cal_a_bp,
             ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus,
             ox_cal_modelled_age_2_σ_uncertainty_cal_a_min,
             region_name
      ) %>%
      mutate(
        t = if_else(!is.na(ox_cal_modelled_age_cal_a_bp) , 
                    ox_cal_modelled_age_cal_a_bp, age_cal_a_bp)
      ) %>%
      mutate(
        Age_2σ_plus = if_else(!is.na(ox_cal_modelled_age_2_σ_uncertainty_cal_a_min) , 
                              ox_cal_modelled_age_2_σ_uncertainty_cal_a_min, 
                              age_2_σ_uncertainty_cal_a_min)
      ) %>%
      mutate(
        Age_2σ_min = if_else(!is.na(ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus) , 
                             ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus, 
                             age_2_σ_uncertainty_cal_a_plus)
      ) %>%
      mutate(
        Z = if_else(type > 0, Z - Z_2σ_min, Z)
      ) %>%
      mutate(
        Z_BgTect = if_else(type > 0, Z_BgTect - Z_BgTect_2σ_min, Z_BgTect)
      )
  }

# Remove tectonic correction from depth when specified
if (correction == 'withBgTect') {
  obs$Z = obs$Z_BgTect
  obs$Z_2σ_min = obs$Z_BgTect_2σ_min
  obs$Z_2σ_plus = obs$Z_BgTect_2σ_plus
  
  # Bounding surfaces not corrected for tectonic subsidence
  HDEM = HDEM_withBgTect
  LDEM = LDEM_withBgTect
  PDEM = PDEM_withBgTect
  
} else {
  HDEM = HDEM_BgTectRem
  LDEM = LDEM_BgTectRem
  PDEM = PDEM_BgTectRem
}

# Ensure that there are no duplicate sets of x-y coordinates, since this will 
# give an error during the Kriging step.
obs <- remove_duplicates(obs)

# Calculate p (normalized age) for each observation.
obs$p = 0
for (i in 1:length(obs$t)){
  if (is.na((A0+obs$t[i])/(A0-A1))) {
    print(paste(obs$ID[i]))
  } else if ((A0+obs$t[i])/(A0-A1) > 1) {
    obs$p[i] = 1
  } else if ((A0+obs$t[i])/(A0-A1) < 0){
    obs$p[i] = 0
  } else {
    obs$p[i] = (A0+obs$t[i])/(A0-A1)}
}

# Process DEMs
PDEM_raster <- create_raster(PDEM, RDnew)
ZA1_raster <- create_raster(HDEM, RDnew)
ZA0_raster <- create_raster(LDEM, RDnew)
Dxy_raster = ZA1_raster - ZA0_raster
q_raster = Dxy_raster/Dmax
VLM_raster <- create_raster(VLM_rate, RDnew)

# Extract point data from raster -----------------------------------------
# Add coordinate system information to observation point dataset
obs_points <- st_as_sf(obs, coords = c("x","y"), crs = RDnew)

# Extract LDEM (=ZA0), HDEM (=ZA1), Dxy and q at observation locations
obs_extract = data.frame(ZA0 = st_extract(ZA0_raster, obs_points),
                         ZA1 = st_extract(ZA1_raster, obs_points),
                         Dxy = st_extract(Dxy_raster, obs_points),
                         q = st_extract(q_raster, obs_points))

# At extracted DEM data to observations
obs$ZA0 = obs_extract[[1]]
obs$ZA1 = obs_extract[[3]]
obs$Dxy = obs_extract[[5]]
obs$q = obs_extract[[7]]

# Calculate relative depth (Zno) of each observation with respect to ZA0 and ZA1
obs$Zno = (obs$Z-obs$ZA0)/obs$Dxy

# Zno between 0 and 1
for (j in 1:length(obs$Z)){
  if (is.na((obs$Z[j]-obs$ZA0[j])/obs$Dxy[j]) == TRUE) {
  } else if(((obs$Z[j]-obs$ZA0[j])/obs$Dxy[j]) > 1) {
    obs$Zno[j] = 1
  } else if (((obs$Z[j]-obs$ZA0[j])/obs$Dxy[j]) < 0){
    obs$Zno[j] = 0
  } else {
    obs$Zno[j] = ((obs$Z[j]-obs$ZA0[j])/obs$Dxy[j])}
}

# Transformation between 0.5 and 1.5 to make sure normalized coordinates are not negative
obs$x_n = 1 -  (mean(grid2D$x) -  obs$x)/(xmax-xmin)
obs$y_n = 1 -  (mean(grid2D$y) -  obs$y)/(ymax-ymin)
obs$age_n = 1 -  (mean(grid3D$age) -  obs$t)/(tmax-tmin)


# Remove raster files that are not being used in the rest of the script
rm(ZA0_raster, ZA1_raster, Dxy_raster, q_raster, PDEM_raster, VLM_raster)


############################# Calculate Zn ################################
# Calculate groundwater trend for observations and for entire grid
index_b <- which(names(params) == "b")

# Calculate parameters
if (coef_lim == "no_limits") { 
  
  a_obs <- calc_a(params, obs$x_n, obs$y_n)
  c_obs <- calc_c(params, obs$x_n, obs$y_n)
  print("No limits used for a and c coefficients")
  name_add_on <- "normal_nolim"
  
} else {
  
  a_obs <- calc_a(params, obs$x_n, obs$y_n)
  for (cell_nr in 1:length(a_obs)) {
    if(a_obs[cell_nr] < 3) {
      a_obs[cell_nr] <- 3
    }
  }
  
  c_obs <- calc_c(params, obs$x_n, obs$y_n)
  for (cell_nr in 1:length(c_obs)) {
    if(c_obs[cell_nr] < 0) {
      c_obs[cell_nr] <- 0
    }
  }
  print("Coefficient a limited to minimal value of 3 and coefficient c limited to minimal value of 0")
  name_add_on <- "minval-a-3_minval-c-0"
}

b <- params[index_b]

# Calculate Zn with the fitted parameters at each observation location
Zn_obs <- sigmoid_function(a_obs, b, c_obs, obs$q, obs$p, obs$x_n, obs$y_n)


# Calculate residuals and transform back to Z(xyt) ------------------------
# Calculate residuals of observations from Z(xyt)_trend and Z(xyt)_obs
obs$Zn_trend = Zn_obs
obs$Zn_res = obs$Zno - obs$Zn_trend

# Transform Zn_trend back to Z(xyt)
obs$Z_trend = obs$ZA0 + obs$Zn_trend * obs$Dxy

# Calculate residuals from Z_trend and Z_obs
obs$Z_res = obs$Z - obs$Z_trend
Z_xyt <- data.frame(Z_trend = obs$Z_trend, age = obs$t,
                    x = obs$x, y = obs$y, Z_res = obs$Z_res)

# Standard deviation of residuals
sigma_Z_res = sd(Z_xyt$Z_res)




########################## Kriging of residuals ###########################
# In this section the residuals that remain after calculating the indicated
# groundwater level as predicted by the trend function are used to determine
# the kriging variogram

# Set coordinate grid for kriging data, using normalized coordinates, 
# since x, y and t are different dimensions!
# This assures that the three axis have the same (normalized) dimension

dx = 1/((xmax-xmin)/xstep)
dy = 1/((ymax-ymin)/ystep)
dt = 1/((tmin-tmax)/tstep)

# Create new dataframe for kriging without NA
Zkrige <- Z_xyt 

Zkrige$x_n = 1 -  (mean(grid2D$x) -  Zkrige$x)/(xmax-xmin)
Zkrige$y_n = 1 -  (mean(grid2D$y) -  Zkrige$y)/(ymax-ymin)
Zkrige$t_n = 1 -  (mean(grid3D$age) -  Zkrige$age)/(tmax-tmin)

# Remove NA fields before kriging
Zkrige <- na.omit(Zkrige) 

# Transfer data to spatial data frame with transformed age-, x-, and y-axes.
coordinates(Zkrige) = ~x_n+y_n+t_n



####################### Kriging variogram fitting #########################
nugget <- ((mean(obs$Z_2σ_min)+mean(obs$Z_2σ_plus))/2)^2
sigma_Z_res <- sd(Zkrige$Z_res)
sill <- sigma_Z_res^2

varmodel_anis <- vgm(psill = sill - nugget, model = "Sph", nugget = nugget, range = .2, 
                     anis = c(90, 90, 90, 0.8, 0.5))


##########
sv_z_age = variogram(Z_res~1, data = Zkrige, cutoff = 0.6, width = 0.05, alpha=c(0), tol.hor = c(180),
                     beta = c(90), tol.ver = c(10))


sv_z_x = variogram(Z_res~1, data = Zkrige, cutoff=0.4, width = 0.03, alpha=c(90), tol.hor=c(20),
                   beta=c(0), tol.ver = c(180))

sv_z_y = variogram(Z_res~1, data = Zkrige, cutoff=0.3, width = 0.021, alpha=c(0),tol.hor=c(40),
                   beta=c(0), tol.ver=c(180))


age_vgram <- plot(sv_z_age, varmodel_anis, plot.numbers = TRUE, main = "Age variogram")
x_vgram <- plot(sv_z_x, varmodel_anis, plot.numbers = TRUE, main = "X variogram")
y_vgram <- plot(sv_z_y, varmodel_anis, plot.numbers = TRUE, main = "Y variogram")

grid.arrange(age_vgram, x_vgram, y_vgram, nrow = 2)


################## Calculate ranges in actual units #######################
range = 0.2
anis_x = 0.8
anis_y = 0.5


range_x <- range * anis_y * (xmax-xmin)/1000 # in [m]
range_y <-range * anis_y * (ymax-ymin)/1000 # in [m]
range_age <-range * -(tmax-tmin) # in [years]

print(range_x)
print(range_y)
print(range_age)

# Transfer output back to normal coordinates
mean(grid2D$x)- (1 - (0.4*0.6))*(xmax-xmin)
mean(grid2D$y)- (1 - Zkrige3D_df$y_n)*(ymax-ymin)
mean(grid3D$age) - (1 - Zkrige3D_df$age_n)*(tmax-tmin)