# Kim de Wit
# The code  presented in this script is part of the 
# “Living on Soft Soils: Subsidence and Society” project 
# (grant no. NWA.1160.18.259), WP1.3

###########################################################################
# Trend fit script paleo-groundwater data Netherlands
###########################################################################
# This script is designed to fit a spatial-temporal trend to Holocene 
# paleo-groundwater level data from the Netherlands. The script uses the 
# the HOLSEA-NL data set as input. The trend that is fitted is designed to 
# have a sigmoid curve shape, which is based on the shape of the relative 
# sea-level rise.

# The input data file needs to be a table with x, y and z coordinates 
# with RD new coordinate projection. 

# NLS fitting method 
# The default algorithm uses a Gauss-Newton algorithm for solving the Least squares

###########################################################################

##### Clear out environment before running the script ############
rm(list = ls())
gc()

# Set working directory and library path ----------------------------------
# Store output in organised folders and files using current month+year (folder)
# and the current day+month+year (output-files)
folder <- format(Sys.time(), "%B_%Y")

if (file.exists(file.path("./results/output", folder, "Trend_fit_output", paste0(Sys.Date())))) {
  cat("The folder already exists")
} else if (file.exists(file.path("./results/output/", folder))) { 
  dir.create(file.path("./results/output/", folder, "Trend_fit_output", paste0(Sys.Date())))
} else {
  dir.create(file.path("./results/output/", folder, "Trend_fit_output", paste0(Sys.Date())), recursive = TRUE)
}

output_loc = paste("./results/output",folder, "Trend_fit_output",
                   paste0(Sys.Date()),sep = "/")


# Loading required packages -----------------------------------------------
# TODO check for dependencies packages
# TODO remove todor package 
library(todor)

library(readxl)
library(snakecase)
library(sp)
library(stars)
library(sf)
library(dplyr)
library(tidyverse)

library(ggplot2)
library(plotly)
library(viridis)
library(ggspatial)
library(grid)
library(gridExtra)
library(metR) # For adding contours

# required for adjusting the headers of the input dataset
source("./src/general_functions.R")

#### Prepare input data ###################################
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

# Select regions to include (or leave out) in trend fit: 
# For the HOLSEA-NL dataset the same regions as in De Wit & Cohen (2024) 
# are used:
# "Noord-Holland", "Rhine-Meuse delta inland", "Flevoland", "Zuid-Holland", 
# "Zeeland", "Waddenzee West", "Waddenzee East"
# Standard all regions are used
regions = c("Noord-Holland", "Rhine-Meuse delta inland", 
            "Flevoland", "Zuid-Holland", "Zeeland",
            "Waddenzee West", "Waddenzee East")


### Define trend fit equation
# Options: "Simplefit","Linearfit_ac","Linearfit_a","Linearfit_c",
# "Quadraticfit_ac","Quadraticfit_a","Quadraticfit_c",
# "Quadraticfit_a_clin", "Quadraticfit_c_alin"



# # a & c linear --------------------------------------------------
# # (b is constant, a and c are linearly dependent on x and y)
# fitmodel_name = "Linearfit_ac"
# 
# startparams = list(a0 = 5, a1 = 5, a2 = 5, 
#                    b = 3 , 
#                    c0 = 0.05, c1 = 0.05, c2 = 0.05) 
# 
# # Define the model formula 
# trend_equation <- Z ~ ZA0 + ((1 - (c0 + c1 * x_n + c2 * y_n)) * 
#                                (1 - exp(-(a0 + a1 * x_n + a2 * y_n) * q * p_fit^b)) + 
#                                (c0 + c1 * x_n + c2 * y_n) * p_fit) * Dxy
# 
# calc_a = function(params, x_n, y_n){
#   params[1] + params[2]*x_n + params[3]*y_n
# }
# 
# calc_c = function(params, x_n, y_n){
#   params[5] + params[6]*x_n + params[7]*y_n
# }



# TODO add trend options to separate config .R file

# a & c linear --------------------------------------------------
# (b is constant, a and c are linearly dependent on x and y)
fitmodel_name = "Linearfit_ac"

startparams = list(a0 = 5, a1 = 5, a2 = 5, 
                   b = 3 , 
                   c0 = 0.05, c1 = 0.05, c2 = 0.05) 

# Define the model formula 
trend_equation <- Z ~ ZA0 + ((1 - (c0 + c1 * x_n + c2 * y_n)) * 
                               (1 - exp(-(a0 + a1 * x_n + a2 * y_n) * q * p_fit^b)) + 
                               (c0 + c1 * x_n + c2 * y_n) * p_fit) * Dxy

calc_a = function(params, x_n, y_n){
  params[1] + params[2]*x_n + params[3]*y_n
}

calc_c = function(params, x_n, y_n){
  params[5] + params[6]*x_n + params[7]*y_n
}




# Create file directory based on input
output_file_loc <- paste(output_loc, region_name, indicator, 
                         correction, fitmodel_name, sep = "/")

if (file.exists(file.path(output_file_loc))) {
  cat("The folder already exists")
} else {
  dir.create(file.path(output_file_loc), recursive = TRUE)
}



# Define Constants -----------------------------------------------------------
Dmean = 11.5  # mean thickness between ZA0 and ZA1
Dmax = 27.1   # max thickness between ZA0 and ZA1
A0 = -10800   # [cal years Before Present] Start time of studied period 
A1 = -1000    # [cal years BP] Start time of studied period 
tmin = -A0
tmax = -A1
tstep = 200   # [years]

t = seq(from = A0, to = A1, by = tstep)  # Time steps
p = (A0-t)/(A0-A1)                       # Normalized time steps

# Grid range - based on extent DEM - in RD new coordinates
xmin = -6000   # [m]
xmax = 280000  # [m]
ymin = 340000  # [m]
ymax = 645000  # [m]
xstep = 1000   # [m]
ystep = 1000   # [m]

# Grid - based on extent DEMs
grid2D = expand.grid(x = seq (from = xmin, to = xmax, by = xstep), 
                     y = seq(from = ymin, to = ymax, by = ystep))
grid3D = expand.grid(x = seq (from = xmin, to = xmax, by = xstep), 
                     y = seq(from = ymin, to = ymax, by = ystep), 
                     age = seq(from = -A0, to = -A1, by = -tstep))

grid2D_points <- st_as_sf(grid2D, coords = c("x","y"), crs = RDnew)


# Shape study area
StudyArea = "data/processed/StudArea_outline.shp"

# Shape provinces NL
provincies = "data/processed/provincies_NL_2021.shp"



# Define bounding surfaces and Pleistocene surface
PDEM = "data/processed/DEMs/PDEM_extended.tif"
PDEM_BgTectRem  = "data/processed/DEMs/PDEM_extended_BgTectCorrected.tif"
HDEM_withBgTect = "data/processed/DEMs/HDEM_extended_NL_v2.tif"
LDEM_withBgTect = "data/processed/DEMs/LDEM_extended_NL_v2.tif"
HDEM_BgTectRem  = "data/processed/DEMs/HDEM_extended_BgTectCorrected.tif"
LDEM_BgTectRem  = "data/processed/DEMs/LDEM_extended_BgTectCorrected.tif"


# Functions #################################
# Define sigmoid function
sigmoid_function = function(a, b, c, q, p, x, y) {
  (1 - c) * (1 - exp(-a * q * p^b)) + (c * p)
}


# Plotting residual data scatter plots against different variables
plot_residualpoints <- function(dataset, variable, residuals, stdv, title,
                                xlabel, ylabel) {
  
  ggplot() +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2*stdv, ymax = 2*stdv),
              fill = "grey", alpha = 0.2) +
    geom_abline(intercept = 0, slope = 0, colour = "grey" , linetype = "longdash", linewidth = 1)+
    geom_segment(aes(x = variable, y = residuals, xend = variable, yend = 0), alpha = 0.2) +
    geom_point(aes(x = variable, y = residuals, fill = residuals), 
               size = 3, shape = 21, colour = "black", stroke = 0.5) +  # Color mapped here
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +  # Colors to use here
    guides(fill = "none") +
    labs(x = xlabel, 
         y = ylabel) +
    geom_text(aes(x = -Inf, y = Inf,
                  label = paste("n = ", length(variable), 
                                "\nsigma = ", round(stdv,3))), 
              hjust = -0.1, vjust = 2,
              fontface = "bold", size = 6) +
    geom_text(aes(x = -Inf, y = 2*stdv,
                  label = paste("2-sigma range")), 
              hjust = -0.1, vjust = 0,
              fontface = "italic", size = 4) +
    ggtitle(title) +
    theme_bw()
  
}

# Plotting age-depth plots at specified locations
plot_agedepth <- function(age, original, fitted, residuals, stdv, title,
                          xlabel, ylabel) {
  ggplot() +
    geom_segment(aes(x = age, y = original, xend = age, yend = fitted), 
                 alpha = 0.6, linetype = "longdash") +
    geom_point(aes(x = age, y = original),
               size = 4, shape = 15, colour = alpha("grey",0.6)) +  
    geom_point(aes(x = age, y = fitted, fill = residuals),
               size = 4, shape = 21, colour = "black", stroke = 0.5) +  # Color mapped here
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +  # Colors to use here
    guides(fill = "none") +
    labs(x = xlabel,
         y = ylabel) +
    geom_text(aes(x = max(age), y = max(original),
                  label = paste("n = ", length(residuals), 
                                "\nsigma = ", round(stdv,3))), 
              hjust = 1, vjust = 1,
              fontface = "bold", size = 6) +
    ggtitle(title) +
    theme_bw()
  
}

# Define color gradient residuals plot
colors_res <- c("#1065AB","#3A93C3","#8EC4DE","#D1E5F0","#F9F9F9",
                "#FEDBC7","#F6A482","#D75F4C","#b31529")

# Parameter map plot function
plot_param_maps <- function(rasterdata, var_name, contour, obs_data, lowerlim, upperlim) { 
  run_info = paste(region_name, indicator, "data", correction, "\n trend fit:", 
                   fitmodel_name, "\n Coast:", coast, "\n n:", n, sep = " ")
  
  ggplot() +  
    geom_stars(data = rasterdata)+
    scale_fill_gradientn(colours = viridis(10, option = "D"), 
                         limits = c(lowerlim, upperlim),
                         name = var_name,
                         na.value = "grey") +
    geom_point(aes(x = obs_data$x, y = obs_data$y), 
               size = 2, shape = 16, colour = alpha("black",0.3)) + 
    geom_contour2(data = AOI_grid, 
                  aes(x = x, y = y, z = contour,
                      label = after_stat(level)),
                  colour = "white", linetype = "dashed")+
    geom_sf(data = provincies_shape, color=alpha("black",0.6), fill = NA, alpha = 1)+
    geom_sf(data = studyarea_shape, color = alpha("red",0.8), alpha = 0, 
            linetype = "dashed") +
    geom_label(aes(x = xmin + 5000, y = ymax - 90000, label = run_info),  # Add a label at a specific point
               fill = "white", color = "black", 
               fontface = "bold.italic", size = 4, hjust = 0) +  # Text box properties
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    # coord_sf(xlim = c(6000, 270000), ylim = c(355000, 630000)) +
    ggtitle(paste(var_name, "-", fitmodel_name,
                  region_name, correction, sep = " "))
  
}

# Z and Zn output map plot function
plot_Z_maps <- function(rasterdata, var_name, contour, countour_interval, 
                        obs_data, lowerlim, upperlim, title) { 
  run_info = paste(region_name, indicator, "data", correction, "\n trend fit:", 
                   fitmodel_name, "\n Coast:", coast, "\n n:", n, sep = " ")
  
  ggplot() +  
    geom_stars(data = rasterdata)+
    scale_fill_gradientn(colours = viridis(10, option = "D"), 
                         limits = c(lowerlim, upperlim),
                         name = var_name,
                         na.value = "grey") +
    # geom_contour2(data = AOI_grid, 
    #               aes(x = x, y = y, z = contour,
    #                   label = after_stat(level)),
    #               breaks = seq(lowerlim, upperlim, by = countour_interval),
    #               colour = "white", linetype = "dashed")+
    geom_point(aes(x = obs_data$x, y = obs_data$y), 
               size = 2, shape = 16, colour = alpha("black",0.3)) + 
    geom_sf(data = provincies_shape, color = alpha("black",0.6), fill = NA, alpha = 1) +
    geom_sf(data = studyarea_shape, color = alpha("red",0.8), alpha = 0, 
            linetype = "dashed") +
    geom_label(aes(x = mean(c(xmax, xmin)), y = 380000, label = run_info),  # Add a label at a specific point
               fill = "white", color = "black", 
               fontface = "bold.italic", size = 3, hjust = 0) +  # Text box properties
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      axis.title = element_blank(),    # Remove all axis labels
      axis.text = element_blank(),     # Remove axis text (tick labels)
      axis.ticks = element_blank()     # Remove axis ticks
    ) + 
    # coord_sf(xlim = c(6000, 270000), ylim = c(355000, 630000)) +
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    ggtitle(title) 

  
}

# Function for creating png file locations
save_as_png <- function(fileloc_base, region_name, indicator, correction,
                        fitmodel_name, variable){ 
  paste(fileloc_base, region_name, indicator, correction, fitmodel_name,
        paste0(fitmodel_name, region_name, "_", variable, ".png"),sep = "/")
}

# Function for saving png files
save_map_pngs <-  function(filename, width, height, plot_name) {
  png(file = filename, width = width, height = height)
  print({
    plot_name +   
      annotation_scale(location = "br")+
      annotation_north_arrow(location = "br", which_north = "true",
                             pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
                             style = north_arrow_nautical()
      ) +
      theme(legend.title = element_text(colour = "black", size = 21, face = "bold"),
            legend.text = element_text(colour = "black", size = 20),
            legend.position = c(.02, .98),
            legend.justification = c("left", "top"),
            legend.box.just = "left",
            legend.margin = margin(2, 6, 8, 6),
            legend.background = element_rect(fill = alpha("white", 0.8),
                                             linewidth = 0.5, linetype = "solid", 
                                             colour ="black")
      )+
      theme(plot.title = element_text(hjust = 0, size = 26, face = "bold"),
            axis.title = element_text(size = 22,face = "bold"),
            axis.text = element_text(size = 21,face ="bold"))
  })
  
  dev.off()
  
}

##### file location for results trend fit
pdf_filename <- paste(output_loc, region_name, indicator, correction, fitmodel_name, 
                      paste0(indicator, correction, fitmodel_name, "_", region_name, paste0(Sys.Date()), 
                             ".pdf"),sep = "/")
txt_filename <- paste(output_loc, region_name, indicator, correction, fitmodel_name,
                      paste0(indicator, correction, fitmodel_name, "_", region_name, paste0(Sys.Date()),
                             "_info.txt"),sep = "/")


# Process input data ---------------------------------------------------------
# Select observation data from HOLSEA data set based on settings
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
                              ox_cal_modelled_age_2_σ_uncertainty_cal_a_min, age_2_σ_uncertainty_cal_a_min)
      ) %>%
      mutate(
        Age_2σ_min = if_else(!is.na(ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus) , 
                             ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus, age_2_σ_uncertainty_cal_a_plus)
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
                              ox_cal_modelled_age_2_σ_uncertainty_cal_a_min, age_2_σ_uncertainty_cal_a_min)
      ) %>%
      mutate(
        Age_2σ_min = if_else(!is.na(ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus) , 
                             ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus, age_2_σ_uncertainty_cal_a_plus)
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
  
  # Bounding surfaces not corrected for tectono-sedimentary subsidence
  HDEM = HDEM_withBgTect
  LDEM = LDEM_withBgTect
  
} else {
  HDEM = HDEM_BgTectRem
  LDEM = LDEM_BgTectRem
}


# Calculate p for each observation. 
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


# Load supporting input files
studyarea_shape = st_read(StudyArea)
provincies_shape = st_read(provincies)

PDEM_raster <- create_raster(PDEM, RDnew)
ZA1_raster <- create_raster(HDEM_withBgTect, RDnew)
ZA0_raster <- create_raster(LDEM_withBgTect, RDnew)
Dxy_raster = ZA1_raster - ZA0_raster
q_raster = Dxy_raster/Dmax

# Extract point data from raster -----------------------------------------
# Add coordinate system information to observation point dataset
obs_points <- st_as_sf(obs, coords = c("x","y"), crs = RDnew)

# Extract GWL1000, GWL10800, Dxy and q at observation locations
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

# Extract GWL1000, GWL10800, Dxy and q at grid locations
AOI_grid = data.frame(ZA0  = st_extract(ZA0_raster, grid2D_points),
                      ZA1  = st_extract(ZA1_raster, grid2D_points),
                      Dxy  = st_extract(Dxy_raster, grid2D_points),
                      q    = st_extract(q_raster, grid2D_points),
                      PDEM = st_extract(PDEM_raster, grid2D_points))

# Remove geometry columns and rename to keep logical column names
AOI_grid <- data.frame(ZA0 = AOI_grid[[1]],
                       ZA1 = AOI_grid[[3]],
                       Dxy = AOI_grid[[5]],
                       q = AOI_grid[[7]],
                       PDEM = AOI_grid[[9]])

# Transformation between 0.5 and 1.5 to make sure normalized coordinates are not negative
AOI_grid$x <- grid2D$x
AOI_grid$y <- grid2D$y
AOI_grid$x_n = 1 -  (mean(grid2D$x) -  grid2D$x)/(xmax-xmin)
AOI_grid$y_n = 1 -  (mean(grid2D$y) -  grid2D$y)/(ymax-ymin)

# Select data before trend fitting for given regions and filter for coastal area when requested  
trendData <- obs %>%
  filter(region_name %in% regions) %>%
  mutate(p_fit = p) %>%
  { 
    if (coast == TRUE) {
      filter(., ZA1 < 1)
    } else {
      .  
    }
  }  


######### TREND FIT ####################################################
# # Fitting to the Z (depth in m) value (and indirectly Zn0)
fitmodel <- nls(trend_equation,
                data = trendData,
                start = startparams,
                control = list(maxiter = 500))

# Show model summary and assign modeled parameters
summary(fitmodel)
params = coef(fitmodel)
index_b <- which(names(startparams) == "b")

# Calculate parameters
a <- calc_a(params, trendData$x_n, trendData$y_n)
c <- calc_c(params, trendData$x_n, trendData$y_n)
b <- params[index_b]

# Calculate predicted Zn based on fitted parameters
Zn_fit <- sigmoid_function(a, b, c, trendData$q, trendData$p_fit, 
                           trendData$x_n, trendData$y_n)

# Predicted values
Z_pred <- predict(fitmodel)

#### STATISTICS ############################################################
# Number of observations and parameters
n <- length(trendData$p_fit)
n_params <- length(coef(fitmodel))

# Get the residuals from the model
residuals <- residuals(fitmodel)

# Calculate sigma2 (variance) and the residual standard error
# The variance and standard error (sigma) calculations are adjusted for the 
# number of parameters in the model (degrees of freedom)
sigma2 <- sum(residuals^2) / (n - n_params)
sigma_rse <- sigma(fitmodel)

# Calculating Cook's distance
# Initialize a vector to store Cook's distance
cooks_dist <- numeric(n)

# Predict Z values using the full model (Z_full)
Z_full <- predict(fitmodel, newdata = trendData)

# Initialize a vector to store Z_without_i values
Z_without <- numeric(nrow(trendData))

# Calculate Cook's distance for each observation
# NOTE: Calculation of Cook's distance can give an error. In those cases, leave
# this section out
for (i in 1:n) {

  # Create a data frame excluding the i-th observation
  data_subset <- data.frame(
    Zno = trendData$Zno[-i],
    x_n = trendData$x_n[-i],
    y_n = trendData$y_n[-i],
    Dxy = trendData$Dxy[-i],
    ZA0 = trendData$ZA0[-i],
    Z = trendData$Z[-i],
    q = trendData$q[-i],
    p_fit = trendData$p[-i]
  )

  # Fit the model without the i-th observation
  nls_model_i <- nls(trend_equation,
                     data = data_subset,
                     start = startparams,
                     control = list(maxiter = 500, minFactor = 1e-13))

  # Predict y value for the i-th observation using the fitted model without it
  Z_without_i <- predict(nls_model_i, newdata = trendData[i, , drop = FALSE])

  # Store the prediction in y_without
  Z_without[i] <- Z_without_i

  # Compute Cook's distance
  cooks_dist[i] <- ((Z_full[i] - Z_without_i)^2) / (n_params * sigma2)
}


# Calculate statistical model metrics ----------------------------------------
# Total Sum of Squares (TSS)
TSS <- sum((trendData$Z - mean(trendData$Z))^2)

# Residual Sum of Squares (RSS): This is calculated as the sum of the squared residuals, which is essentially 
RSS <- sum(residuals^2)

# R-squared
R_squared <- 1 - (RSS / TSS)

# Adjusted R-squared
adjusted_R_squared <- 1 - (1 - R_squared) * ((n - 1) / (n - n_params))

# Root Mean Squared Error (RMSE)
RMSE = sqrt(mean(residuals^2))

# Metrics for Zn
# Total Sum of Squares (TSS)
TSS_Zn <- sum((trendData$Zno - mean(trendData$Zno))^2)

# Residual Sum of Squares (RSS)
residualsZn = trendData$Zno - Zn_fit
RSS_Zn <- sum(residualsZn^2)

# R-squared
R_squaredZn <- 1 - (RSS_Zn / TSS_Zn)

# Adjusted R-squared
adjusted_R_squaredZn <- 1 - (1 - R_squaredZn) * ((n - 1) / (n - n_params))

# Root Mean Squared Error (RMSE)
RMSE_Zn = sqrt(mean(residualsZn^2))

# Standard deviation Zn and Z. For sigma_Z the residual standard error is used
# which is adjusted for the degrees of freedom of the model
sigma_Z = sigma_rse
sigma_Zn = sd(residualsZn)


# Add predicted data and residuals to trendData frame for visualization
trendData$Zn <- Zn_fit
trendData$resZn <- residualsZn
trendData$Zpred <- Z_pred
trendData$resZ <- residuals


# Plot Cook's distance -----------------
# Define the threshold for Cook's distance
threshold <- 8/(n - 2 * (n_params + 1))

# Identify the points where Cook's distance exceeds the threshold
exceeding_threshold <- cooks_dist > threshold

# Scale-Location plot input
# Calculate standardized residuals
standardized_residuals <- residuals / sigma_rse

# Calculate the square root of the absolute standardized residuals
sqrt_abs_residuals <- sqrt(abs(standardized_residuals))

# Calculate leverage
# For non-linear models, leverage can be calculated using the hat values.
# Here we approximate this using the hatvalues() function from the linear model (lm).

# Fit an auxiliary linear model to estimate leverage values
auxiliary_model <- lm(Z ~ x_n + y_n, data = trendData)

# Calculate leverage
leverage_values <- hatvalues(auxiliary_model)


######### SAVE OUTPUT ##################################################

# Save model information and statistical information in text file ------------
sink(file = txt_filename)

print(summary(fitmodel)) # Model that was used for trendfitting

if (coast == TRUE) {
  print("Coastal data selected, below HDEM = 1 m")
} else {
  print("Not filtered for coastal data")
}

print(obs_type)

cat(" R squared Z = ", R_squared,"\n",
              "Adjusted R squared Z = ", adjusted_R_squared,"\n",
              "Sigma Z = ", sigma_Z,"\n",
              "RMSE = ", RMSE,"\n",
              "RSS = ", RSS,"\n",
              "R squared Zn = ", R_squaredZn,"\n",
              "Adjusted R squared Zn = ", adjusted_R_squared,"\n",
              "Sigma Zn = ", sigma_Zn,"\n",
              "RMSE Zn = ", RMSE_Zn,"\n",
              "RSS Zn = ", RSS_Zn,"\n",
              "n = ", length(trendData$x),
              "\n", "correction = ", correction)

sink()

graphics.off()

# Save observation data as csv 
write.csv(obs, paste(output_loc, region_name, indicator, correction, fitmodel_name,
                     "/obs_data.csv", sep = "/"))

write.csv(trendData, paste(output_loc, region_name, indicator, correction, fitmodel_name,
                     "/trendData.csv", sep = "/"))

# Save results parameter fitting as PDF, pngs and text files ----------------
# Open pdf file
pdf(file = pdf_filename)

layout_matrix_2 <- matrix(1:4, ncol = 2)
layout(layout_matrix_2)
# Standard figures of input data
plot(obs$t, obs$Z, 
     xlab = "Age [cal. yrs BP]", 
     ylab = "Depth Z [m]",
     main = "Observations depth with age")
plot(obs$t, obs$Dxy, 
     xlab = "Age [cal. yrs BP]", 
     ylab = "Tickness Holocene wedge D [m]",
     main = "Thickness Holocene wedge obs with age")
plot(obs$t, obs$Zno, 
     xlab = "Age [cal. yrs BP]", 
     ylab = "Normalized depth Zno [m]",
     main = "Normalized observation depth with age")
plot(obs$Z, obs$Zno, 
     xlab = "Depth Z [m]", 
     ylab = "Normalized depth Zno [m]",
     main = "Actual observation depth vs normalized depth")

grid.arrange(
  ggplot(data = obs)+
  geom_point(aes(x = t, y = Zno, color = x))+
  scale_color_gradientn(colours = viridis(256, option = "D"),
                        name = "x-coordinate") + 
  ggtitle("Spread of observations East-West"),

ggplot(data = obs)+
  geom_point(aes(x = t, y = Zno, color = y))+
  scale_color_gradientn(colours = viridis(256, option = "D"),
                        name = "y-coordinate") + 
  ggtitle("Spread of observations North-South"),
ncol = 2
)

# Plot diagnostic plots in one overview
layout_matrix_1 <- matrix(1:6, ncol = 2)
layout(layout_matrix_1)

plot(trendData$Z, Z_pred, pch = 18, col = alpha("black",0.6),
     xlab = "Index point depth [m]", 
     ylab = "Fitted depth [m]",
     main = "Fitted vs Observed")
abline(0, 1, lty = 2)
abline(0 + sigma_Z, 1, lty = 3, col = alpha("black",0.6))
abline(0 - sigma_Z, 1, lty = 3, col = alpha("black",0.6))

# Q-Q plot to check if residuals are normal distributed
qqnorm(residuals, pch = 1)
qqline(residuals, col = "steelblue", lwd = 2)

# Create the Cook's distance plot
plot(cooks_dist, type = "h", lwd = 2,
     xlab = "Observation Index", 
     ylab = "Cook's Distance",
     main = "Cook's Distance Plot")

# Add a horizontal line for the threshold
abline(h = threshold, col = "blue", lty = 2)
abline(h = 4/(n - n_params), col = "red", lty = 2)  # Threshold line


# Add labels only for points that exceed the threshold
if (length(which(exceeding_threshold)) == 0) {
  text(x = n - 100,
       y = max(cooks_dist),
       label = "All points below treshold")
} else {
text(x = which(exceeding_threshold), 
     y = cooks_dist[exceeding_threshold], 
     labels = trendData$Name[which(exceeding_threshold)], 
     pos = 4, cex = 0.6)
}

# Residuals vs Fitted plot
plot(Z_pred, residuals, 
     xlab = "Fitted values", 
     ylab = "Residuals",
     main = "Residuals vs Fitted")
lines(lowess(Z_pred, residuals), col = "red", lwd = 2)

# Scale-Location plot
plot(Z_pred, sqrt_abs_residuals, 
     xlab = "Fitted Values", 
     ylab = "√|Standardized Residuals|",
     main = "Scale-Location Plot")

# Add a smooth line to help identify patterns
abline(h = 0, col = "blue", lty = 3)  # Optional: horizontal reference line
lines(lowess(Z_pred, sqrt_abs_residuals), col = "red", lwd = 2)

# Residuals vs. Leverage plot
plot(leverage_values, standardized_residuals, 
     xlab = "Leverage", 
     ylab = "Standardized Residuals",
     main = "Standardized Residuals vs. Leverage")

# 4. Add a horizontal line at 0 to help identify patterns
abline(h = 0, col = "lightblue", lty = 3)
lines(lowess(leverage_values, standardized_residuals), col = "red", lwd = 2)

if (length(which(exceeding_threshold)) == 0) {
  text(x = n - 100,
       y = max(cooks_dist),
       label = "All points below treshold")
} else {
  # Add labels only for points that exceed the threshold
  text(x = leverage_values[exceeding_threshold], 
       y = standardized_residuals[exceeding_threshold], 
       labels = trendData$Name[which(exceeding_threshold)], 
       pos = 4, cex = 0.6)
}


layout_matrix_0 <- matrix(1)
layout(layout_matrix_0)

# Draw plots
## Residuals Zn
# Histogram plot
ggplot()+
  geom_histogram(aes(y = after_stat(density), x = residualsZn), 
                 alpha = 0.6, binwidth = 0.05, 
                 fill =  "#69b3a2", 
                 color = "#69b3a2", 
                 boundary = 0) +
  stat_function(fun = dnorm, args = list(mean = mean(residualsZn), 
                                         sd = sd(residualsZn)), 
                color = "#404080", linewidth = 1) +
  xlab("Relative vertical deviation Zn") +
  ggtitle("Residuals distribution of Zn") +
  theme_bw()


# Scatterplots
plot_residualpoints(trendData, Zn_fit, residualsZn, sigma_Zn,
                    "Distribution residuals with Zn", "Zn", 
                    "Relative vertical deviation Zn")

plot_residualpoints(trendData, trendData$x, residualsZn, sigma_Zn,
                    "Distribution residuals with x coordinate", "x [m]", 
                    "Relative vertical deviation Zn")

plot_residualpoints(trendData, trendData$y, residualsZn, sigma_Zn,
                    "Distribution residuals with y coordinate", "y [m]", 
                    "Relative vertical deviation Zn")

plot_residualpoints(trendData, trendData$Dxy, residualsZn, sigma_Zn,
                    "Distribution residuals with thickness", "Dxy [m]", 
                    "Relative vertical deviation Zn")

plot_residualpoints(trendData, trendData$t, residualsZn, sigma_Zn,
                    "Distribution residuals with age", "Age [years cal BP]", 
                    "Relative vertical deviation Zn")

plot_agedepth(trendData$t, trendData$Zno, Zn_fit, trendData$resZn, sigma_Zn, 
              "Age-depth plot actual and fitted Zn", "Age", "Zn")

ggplot() +  
  geom_sf(data = provincies_shape, color = alpha("black",0.6), fill = "grey88", alpha = 1) +
  geom_point(aes(x = trendData$x, y = trendData$y, fill = trendData$resZn, 
                 size = abs(trendData$resZn)), 
             shape = 21, colour = "black", stroke = 0.5) +  # Color mapped here
  scale_fill_gradientn(colors = colors_res,
                       # limits = c(-4, 5), 
                       name = "Residuals") +
  guides(size = "none") +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.title = element_blank(),    # Remove all axis labels
    axis.text = element_blank(),     # Remove axis text (tick labels)
    axis.ticks = element_blank(),     # Remove axis ticks
    legend.position = "bottom"
  ) + 
  ggtitle("Map of Zn residuals") +
  coord_sf(xlim = c(6000, 270000), ylim = c(355000, 630000))


# Histogram plot
ggplot()+
  geom_histogram(aes(y = after_stat(density), x = residuals), 
                 alpha = 0.6, binwidth = 0.5, 
                 fill =  "#69b3a2", 
                 color = "#69b3a2", 
                 boundary = 0) +
  stat_function(fun = dnorm, args = list(mean = mean(residuals), 
                                         sd = sd(residuals)), 
                color = "#404080", linewidth = 1) +
  xlab("Vertical deviation Z [m]") +
  ggtitle("Residuals distribution of Z") +
  theme_bw()

# Scatter plots
plot_residualpoints(trendData, trendData$Z, residuals, sigma_Z,
                    "Distribution residuals with Z", "Z [m]", 
                    "Vertical deviation Z [m]")

plot_residualpoints(trendData, trendData$x, residuals, sigma_Z,
                    "Distribution residuals with x coordinate", "x [m]", 
                    "Vertical deviation Z [m]")

plot_residualpoints(trendData, trendData$y, residuals, sigma_Z,
                    "Distribution residuals with y coordinate", "y [m]", 
                    "Vertical deviation Z [m]")

plot_residualpoints(trendData, trendData$Dxy, residuals, sigma_Z,
                    "Distribution residuals with thickness", "Dxy [m]", 
                    "Vertical deviation Z [m]")

plot_residualpoints(trendData, trendData$t, residuals, sigma_Z,
                    "Distribution residuals with age", "Age [years cal BP]", 
                    "Vertical deviation Z [m]")


plot_agedepth(trendData$t, trendData$Z, trendData$Zpred, trendData$resZ, sigma_Z, 
              "Age-depth plot actual and fitted Z",
              "Age", "Z [m]")

ggplot() +  
  geom_sf(data = provincies_shape, color = alpha("black",0.6), fill = "grey88", alpha = 1) +
  geom_point(aes(x = trendData$x, y = trendData$y, fill = trendData$resZ, 
                 size = abs(trendData$resZ)), 
             # size = 5, 
             shape = 21, colour = "black", stroke = 0.5) +  # Color mapped here
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red") +  # Colors to use here
  scale_fill_gradientn(colors = colors_res,
                       # limits = c(-4, 5), 
                       name = "Residuals [m]") +
  # breaks = seq(-4, 4, 1)) +  # Colors to use here
  scale_size(breaks = seq(0, 4, 0.5)) +
  guides(size = "none") +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.title = element_blank(),    # Remove all axis labels
    axis.text = element_blank(),     # Remove axis text (tick labels)
    axis.ticks = element_blank(),     # Remove axis ticks
    legend.position = "bottom"
  ) + 
  ggtitle("Map of Z residuals") +
  coord_sf(xlim = c(6000, 270000), ylim = c(355000, 630000))


plot(trendData$Z, trendData$Zpred, pch = 18, col=alpha("black",0.6))
abline(0, 1, lty = 2)
abline(0 + sigma_Z, 1, lty = 3, col = alpha("black",0.6))
abline(0 - sigma_Z, 1, lty = 3, col = alpha("black",0.6))

plot(trendData$Zno, trendData$Zn, pch = 18, col=alpha("black",0.6))
abline(0, 1, lty = 2)
abline(0 + sigma_Zn, 1, lty = 3, col = alpha("black",0.6))
abline(0 - sigma_Zn, 1, lty = 3, col = alpha("black",0.6))
dev.off()

graphics.off()



# Map plots of spatial patterns parameters --------------------------------
# Calculate parameters a and c separate for entire grid
AOI_grid$a <- calc_a(params, AOI_grid$x_n, AOI_grid$y_n)
AOI_grid$c <- calc_c(params, AOI_grid$x_n, AOI_grid$y_n)
AOI_grid$cDxy <- (AOI_grid$c*AOI_grid$Dxy)/(A0-A1)*1000 # approximation of linear movement rate [m/kyr]
AOI_grid$aq <- (AOI_grid$a*AOI_grid$q)

param_raster <- st_as_sf(AOI_grid, coords = c("x","y"), crs = RDnew)
param_raster = st_rasterize(param_raster, dx = 1000, dy = 1000)

# Standard parameter distribution maps
plot_a <- plot_param_maps(param_raster[8], "Coefficient a", AOI_grid$a, obs, 0, 30)
plot_aq <- plot_param_maps(param_raster[8]*param_raster[4], "a(x,y)*q(x,y)", AOI_grid$aq, obs, 0, 30)
plot_c <- plot_param_maps(param_raster[9], "Coefficient c", AOI_grid$c, obs, -0.5, 0.5)
plot_cDxy <- plot_param_maps(param_raster[10], "[m/kyr]", AOI_grid$cDxy, obs, -0.5, 0.5)

# Specify output location
png_a_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "a")
png_aq_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "aq")
png_c_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "c")
png_cDxy_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "cDxy")

# Save parameter plots
save_map_pngs(png_a_filename, 800, 800, plot_a) 
save_map_pngs(png_aq_filename, 800, 800, plot_aq) 
save_map_pngs(png_c_filename, 800, 800, plot_c) 
save_map_pngs(png_cDxy_filename, 800, 800, plot_cDxy) 

# Map plots of spatial patterns Zn and Z --------------------------------
# Zn with equal thickness of Holocene sediment of 11.75 m
Zn_sameD_df <- data.frame(
  Zn_2k = sigmoid_function(AOI_grid$a, b, AOI_grid$c, 0.5, 
                           p[which(t == -2000)], AOI_grid$x_n, AOI_grid$y_n),
  Zn_4k = sigmoid_function(AOI_grid$a, b, AOI_grid$c, 0.5, 
                           p[which(t == -4000)], AOI_grid$x_n, AOI_grid$y_n),
  Zn_6k = sigmoid_function(AOI_grid$a, b, AOI_grid$c, 0.5, 
                           p[which(t == -6000)], AOI_grid$x_n, AOI_grid$y_n),
  Zn_8k = sigmoid_function(AOI_grid$a, b, AOI_grid$c, 0.5, 
                           p[which(t == -8000)], AOI_grid$x_n, AOI_grid$y_n),
  Zn_10k = sigmoid_function(AOI_grid$a, b, AOI_grid$c, 0.5, 
                            p[which(t == -10000)], AOI_grid$x_n, AOI_grid$y_n),
  x = grid2D$x, 
  y = grid2D$y
  )

# Zn with actual thickness of Holocene sediments
Zn_actualD_df <- data.frame(
  Zn_2k = sigmoid_function(AOI_grid$a, b, AOI_grid$c, AOI_grid$q, 
                           p[which(t == -2000)], AOI_grid$x_n, AOI_grid$y_n),
  Zn_4k = sigmoid_function(AOI_grid$a, b, AOI_grid$c, AOI_grid$q, 
                           p[which(t == -4000)], AOI_grid$x_n, AOI_grid$y_n),
  Zn_6k = sigmoid_function(AOI_grid$a, b, AOI_grid$c, AOI_grid$q, 
                           p[which(t == -6000)], AOI_grid$x_n, AOI_grid$y_n),
  Zn_8k = sigmoid_function(AOI_grid$a, b, AOI_grid$c, AOI_grid$q, 
                           p[which(t == -8000)], AOI_grid$x_n, AOI_grid$y_n),
  Zn_10k = sigmoid_function(AOI_grid$a, b, AOI_grid$c, AOI_grid$q, 
                            p[which(t == -10000)], AOI_grid$x_n, AOI_grid$y_n),
  x = grid2D$x, 
  y = grid2D$y
)

# Z with actual thickness of Holocene sediments
Z_df <- data.frame(
  Z_2K = AOI_grid$ZA0 + sigmoid_function(AOI_grid$a, b, AOI_grid$c, AOI_grid$q, 
                                         p[which(t == -2000)], 
                                         AOI_grid$x_n,  AOI_grid$y_n) * AOI_grid$Dxy,
  Z_4k = AOI_grid$ZA0 + sigmoid_function(AOI_grid$a, b, AOI_grid$c, AOI_grid$q, 
                                         p[which(t == -4000)], 
                                         AOI_grid$x_n,  AOI_grid$y_n) * AOI_grid$Dxy,
  Z_6k = AOI_grid$ZA0 + sigmoid_function(AOI_grid$a, b, AOI_grid$c, AOI_grid$q, 
                                         p[which(t == -6000)], 
                                         AOI_grid$x_n,  AOI_grid$y_n) * AOI_grid$Dxy,
  Z_8k = AOI_grid$ZA0 + sigmoid_function(AOI_grid$a, b, AOI_grid$c, AOI_grid$q, 
                                         p[which(t == -8000)], 
                                         AOI_grid$x_n,  AOI_grid$y_n) * AOI_grid$Dxy,
  Z_10k = AOI_grid$ZA0 + sigmoid_function(AOI_grid$a, b, AOI_grid$c, AOI_grid$q, 
                                          p[which(t == -10000)], 
                                          AOI_grid$x_n,  AOI_grid$y_n) * AOI_grid$Dxy,
  x = grid2D$x, 
  y = grid2D$y
)


# Zn with one thickness
Zn_sameD_sf <- st_as_sf(Zn_sameD_df, coords = c("x","y"), crs = RDnew)
Zn_sameD_raster <- st_rasterize(Zn_sameD_sf, dx = 1000, dy = 1000)

# Zn with actual thickness
Zn_actualD_sf <- st_as_sf(Zn_actualD_df, coords = c("x","y"), crs = RDnew)
Zn_actualD_raster <- st_rasterize(Zn_actualD_sf, dx = 1000, dy = 1000)

# Zn with actual thickness
Z_sf <- st_as_sf(Z_df, coords = c("x","y"), crs = RDnew)
Z_raster <- st_rasterize(Z_sf, dx = 1000, dy = 1000)

# Overview plots of Zn with same D with color bar scale from 0 to 1
plot_ZnsameD_2k <- plot_Z_maps(Zn_sameD_raster[1], "Zn_sameD", Zn_sameD_df$Zn_2k, 0.5, obs, 0, 1, "2000 BP")
plot_ZnsameD_4k <- plot_Z_maps(Zn_sameD_raster[2], "Zn_sameD", Zn_sameD_df$Zn_4k, 0.5, obs, 0, 1, "4000 BP")
plot_ZnsameD_6k <- plot_Z_maps(Zn_sameD_raster[3], "Zn_sameD", Zn_sameD_df$Zn_6k, 0.5, obs, 0, 1, "6000 BP")
plot_ZnsameD_8k <- plot_Z_maps(Zn_sameD_raster[4], "Zn_sameD", Zn_sameD_df$Zn_8k, 0.5, obs, 0, 1, "8000 BP")
plot_ZnsameD_10k <- plot_Z_maps(Zn_sameD_raster[5], "Zn_sameD", Zn_sameD_df$Zn_10k, 0.5, obs, 0, 1, "10000 BP")

# Overview plots of Zn with same D with variable color bar scale
plot_ZnsameD_nolim_2k <- plot_Z_maps(Zn_sameD_raster[1], "Zn_sameD", Zn_sameD_df$Zn_2k, 0.5, obs,
                               min(Zn_sameD_df$Zn_2k), max(Zn_sameD_df$Zn_2k), "2000 BP")
plot_ZnsameD_nolim_4k <- plot_Z_maps(Zn_sameD_raster[2], "Zn_sameD", Zn_sameD_df$Zn_4k, 0.5, obs,
                                     min(Zn_sameD_df$Zn_4k), max(Zn_sameD_df$Zn_4k), "4000 BP")
plot_ZnsameD_nolim_6k <- plot_Z_maps(Zn_sameD_raster[3], "Zn_sameD", Zn_sameD_df$Zn_6k, 0.5, obs,
                                     min(Zn_sameD_df$Zn_6k), max(Zn_sameD_df$Zn_6k), "6000 BP")
plot_ZnsameD_nolim_8k <- plot_Z_maps(Zn_sameD_raster[4], "Zn_sameD", Zn_sameD_df$Zn_8k, 0.5, obs,
                                     min(Zn_sameD_df$Zn_8k), max(Zn_sameD_df$Zn_8k), "8000 BP")
plot_ZnsameD_nolim_10k <- plot_Z_maps(Zn_sameD_raster[5], "Zn_sameD", Zn_sameD_df$Zn_10k, 0.5, obs,
                                      min(Zn_sameD_df$Zn_10k), max(Zn_sameD_df$Zn_10k), "10000 BP")

# Overview plots of Zn with actual D.
plot_ZnactualD_2k <- plot_Z_maps(Zn_actualD_raster[1], "Zn_actualD", Zn_actualD_df$Zn_2k, 0.5, obs, 0, 1, "2000 BP")
plot_ZnactualD_4k <- plot_Z_maps(Zn_actualD_raster[2], "Zn_actualD", Zn_actualD_df$Zn_4k, 0.5, obs, 0, 1, "4000 BP")
plot_ZnactualD_6k <- plot_Z_maps(Zn_actualD_raster[3], "Zn_actualD", Zn_actualD_df$Zn_6k, 0.5, obs, 0, 1, "6000 BP")
plot_ZnactualD_8k <- plot_Z_maps(Zn_actualD_raster[4], "Zn_actualD", Zn_actualD_df$Zn_8k, 0.5, obs, 0, 1, "8000 BP")
plot_ZnactualD_10k <- plot_Z_maps(Zn_actualD_raster[5], "Zn_actualD", Zn_actualD_df$Zn_10k, 0.5, obs, 0, 1, "10000 BP")

# Overview plots of Z with actual D.
plot_Z_2k <- plot_Z_maps(Z_raster[1], "Z [m]", Z_df$Z_2K, 1, obs, -30, 5, "2000 BP")
plot_Z_4k <- plot_Z_maps(Z_raster[2], "Z [m]", Z_df$Z_4k, 1, obs, -30, 5, "4000 BP")
plot_Z_6k <- plot_Z_maps(Z_raster[3], "Z [m]", Z_df$Z_6k, 1, obs, -30, 5, "6000 BP")
plot_Z_8k <- plot_Z_maps(Z_raster[4], "Z [m]", Z_df$Z_8k, 1, obs, -30, 5, "8000 BP")
plot_Z_10k <- plot_Z_maps(Z_raster[5], "Z [m]", Z_df$Z_10k, 1, obs, -30, 5, "10000 BP")

# Specify output location
png_ZnsameDxy_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "ZnsameDxy")
png_ZnsameD_nolim_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "ZnsameD_nolim")
png_Zn_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name,  "Zn_actualD")
png_Z_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name,  "Z_actualD")


# Save Z map plots
save_map_pngs(png_ZnsameDxy_filename, 1600, 800, 
              grid.arrange(plot_ZnsameD_10k, plot_ZnsameD_8k, plot_ZnsameD_6k, 
                           plot_ZnsameD_4k, plot_ZnsameD_2k, nrow = 2)) 

save_map_pngs(png_ZnsameD_nolim_filename, 1600, 800, 
              grid.arrange(plot_ZnsameD_nolim_10k, plot_ZnsameD_nolim_8k, plot_ZnsameD_nolim_6k, 
                           plot_ZnsameD_nolim_4k, plot_ZnsameD_nolim_2k, nrow = 2)) 

save_map_pngs(png_Zn_filename, 1600, 800, 
              grid.arrange(plot_ZnactualD_10k, plot_ZnactualD_8k, plot_ZnactualD_6k, 
                           plot_ZnactualD_4k, plot_ZnactualD_2k, nrow = 2)) 

save_map_pngs(png_Z_filename, 1600, 800, 
              grid.arrange(plot_Z_10k, plot_Z_8k, plot_Z_6k, 
                           plot_Z_4k, plot_Z_2k, nrow = 2)) 



print(pdf_filename)
print(txt_filename)
