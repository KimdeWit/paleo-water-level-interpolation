# Kim de Wit
# The code  presented in this script is part of the 
# “Living on Soft Soils: Subsidence and Society” project 
# (grant no. NWA.1160.18.259), WP1.3

###########################################################################
# Interpolation visualization script
###########################################################################
# Script for visualizing paleo-water level interpolation output, useful for 
# analyzing the output.
# This script should be ran directly after "Mainscript-GWL-interpolation.R", 
# to ensure all required input data is loaded.

################################ Libraries ################################
library(grid)
library(gridExtra)
library(ggplot2)
library(plotly)
library(viridis)
library(ggspatial)
library(metR) # For adding contours


############################## Define Constants ###########################
# Shapes other countries for the background
Belgium <- readRDS("./data/raw/Country_shapefiles/gadm36_BEL_0_sp.rds")
Belgium <- st_as_sf(Belgium, crs = RDnew)
Germany <- readRDS("./data/raw/Country_shapefiles/gadm36_DEU_0_sp.rds")
Germany <- st_as_sf(Germany, crs = RDnew)


# Shape study area
StudyArea = "data/processed/StudArea_outline.shp"

# Shape provinces NL
provincies = "data/processed/provincies_NL_2021.shp"

# Load supporting input files
studyarea_shape = st_read(StudyArea)
provincies_shape = st_read(provincies)


PDEM_raster <- create_raster(PDEM, RDnew)
ZA1_raster <- create_raster(HDEM, RDnew)
ZA0_raster <- create_raster(LDEM, RDnew)
# Dxy_raster = ZA1_raster - ZA0_raster
# q_raster = Dxy_raster/Dmax
# VLM_raster <- create_raster(VLM_rate, RDnew)



############################### Functions #################################
save_as_png <- function(fileloc_base, indicator, correction, 
                        fitmodel_name, variable){ 
  paste(fileloc_base, indicator, correction,
        paste0(fitmodel_name, "_", variable, ".png"),sep = "/")
}

save_ad_plots <-  function(filename, width, height, plot_name) {
  png(file = filename, width = width, height = height)
  print({
    plot_name
  })
  
  dev.off()
  
}

save_map_plots <-  function(filename, width, height, plot_name) {
  png(file = filename, width = width, height = height)
  print({
    plot_name
  })
  
  dev.off()
  
}

plot_interpolation_maps <- function(layer, name, timesteps, obs_data) {
  
  
  map_color <- c("H", "D", "D", "H", "B", "A", "H")
  col_dir <- c(1, 1, 1, -1, 1, 1)
  
  if (layer == 1 | layer == 6) {
    # Zn & variance
    legend_breaks <- c(0, 0.5, 1)
    legend_labels <- c("0", "0.5", "1")
    legend_limits <- c(0, 1)
  } else if (layer == 2) {
    # PDEM
    legend_breaks <- c(-30, -10, 0, 7)
    legend_labels <- c("-30", "-10", "0", "7")
    legend_limits <- c(-30, 7)
  } else if (layer == 3) {
    # VLM_rate
    legend_breaks <- c(-0.2, -0.1, 0)
    legend_labels <- c("-0.2", "-0.1", "0")
    legend_limits <- c(-0.2, 0)
  } else if (layer == 5) {
    # Kriging prediction
    legend_breaks <- c(-5, 0 , 5)
    legend_labels <- c("-5", "0", "5")
    legend_limits <- c(-5, 5)
  } else {
    # Z_Trend, Zxyt
    legend_breaks <- c(-30, -20, -10, 0, 10)
    legend_labels <- c("-30", "-20", "-10", "0", "10")
    legend_limits <- c(-30, 10)
  }
  
  
  timeslice_maps <- ggplot() + 
    geom_sf(data = provincies_shape, color=alpha("black",0.2), 
            fill = "grey", alpha = 1) +
    geom_stars(data = time_extract[layer, , , timesteps], alpha = 1) +
    scale_fill_gradientn(colours = viridis(256, option = map_color[layer],
                                           direction = 1),
                         name = name,
                         breaks = legend_breaks,
                         labels = legend_labels,
                         limits = legend_limits,
                         na.value = NA)+
    # facet_wrap("age", nrow = 2) +
    # facet_wrap("age") +
    ggtitle(paste("interpolated", indicator, "data", correction)) +
    geom_sf(data = provincies_shape, color=alpha("black",0.2), fill = NA, alpha = 0.5) +
    geom_sf(data = Belgium, color=alpha("black",1), fill = "grey50", alpha = 1)+
    geom_sf(data = Germany, color=alpha("black",1), fill = "grey50", alpha = 1)+
    geom_sf(data = provincies_shape, color=alpha("black",0.2), fill = NA, alpha = 0.5) +
    geom_point(data = obs_data, aes(x = x, y = y),
               size = 1, shape = 16, colour = alpha("black",0.4)) +
    facet_wrap("age", nrow = 2) +
    coord_sf(xlim = c(6000, 270000), ylim = c(365000, 630000)) +
    guides(fill = guide_colourbar(title.position = "right"))+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 20, face = "bold")) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(size = 30, face = "bold")) +
    theme(legend.position = "right", legend.key.height = unit(2, "cm"),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 30, angle = -90, face = "bold"),
          legend.title.align = 0.5,
          legend.direction = "vertical",
          strip.text.x = element_text(size = 20, color = "black", face = "bold.italic"),
          plot.title = element_text(hjust = 0,size = 20, face = "bold"))
  
  timeslice_maps
  
  
}


ad_plot_single <-  function(x, y, name) {
  plot_loc <- data.frame(x = x, y = y, location = name)
  
  # Extract GWL1000, GWL10800, Dxy and q at observation locations
  plot_loc <- st_as_sf(plot_loc, coords = c("x","y"), crs = RDnew)
  
  
  #### Adding hlines for DEMs
  plot_loc_points <- st_as_sf(plot_loc, coords = c("x","y"), crs = RDnew)
  
  plot_loc_extract <- data.frame(ZA0 = st_extract(ZA0_raster, plot_loc_points),
                                 ZA1 = st_extract(ZA1_raster, plot_loc_points),
                                 # PDEM = st_extract(star_all[2,,,1], plot_loc_points))
                                 PDEM = st_extract(star_all[2,,,1], plot_loc_points))
  
  
  
  # Extract the data at each specified location from the interpolation dataset
  GWLcurve_extract <- st_extract(star_all, plot_loc)
  
  # Store extracted data, timesteps and location name in new data frame
  GWLcurves_data <- data.frame()
  for (k in 1){
    locdata =   as.data.frame(GWLcurve_extract[,k,])
    locdata$age = -t
    locdata$Location = plot_loc$location[k]
    if (k == 1){
      GWLcurves_data = locdata
    } else {
      GWLcurves_data = rbind(locdata, GWLcurves_data)
    }
  }
  
  
  obs_loc <- obs %>%
    filter(grepl(name, Name))
  
  SLIP_data <- obs_loc %>%
    filter(type %in% "0")
  
  ULD_data <- obs_loc %>%
    filter(type %in% c("1", "2", "3"))
  
  run_info = paste(indicator, "data", correction, "\n trend fit:",
                   fitmodel_name, "\n Sigma:", round(sigma_Z_res, 2), "m", 
                   "\n n:", length(obs$ID), sep = " ")
  
  type_colors <- c("1" = "#0072B2", "2" = "lightseagreen", "3" = "#D55E00")
  
  ggplot() +
    geom_rect(aes(xmin = -Inf, xmax = Inf,   
                  ymin = -Inf, ymax = plot_loc_extract[[5]]),  
              fill = "burlywood", alpha = 0.3    
    ) +
    geom_hline(yintercept = plot_loc_extract[[5]],    # PDEM
               color = "burlywood",  , alpha = 0.5,
               linetype = "solid", 
               linewidth = 1) + 
    geom_text(data = NULL, aes(x = A1, y = plot_loc_extract[[5]], 
                               label = "Pleistocene surface"), 
              size = 6, hjust = 1, vjust = -0.2) +
    geom_hline(yintercept = plot_loc_extract[[1]],    # LDEM
               color = "navyblue", 
               linetype = "dotted", 
               linewidth = 1) + 
    geom_text(data = NULL, aes(x = A1, y = plot_loc_extract[[1]], 
                               label = "Lowstand"), 
              size = 6, hjust = 1, vjust = -0.2) +
    geom_hline(yintercept = plot_loc_extract[[3]],    # HDEM
               color = "navyblue", 
               linetype = "dotted", 
               linewidth = 1) + 
    geom_text(data = NULL, aes(x = A0, y = plot_loc_extract[[3]], 
                               label = "Highstand"), 
              size = 6, hjust = 0, vjust = -0.2) +
    
    geom_ribbon(data = GWLcurves_data, aes(x = -age, ymin = Ztrend - 2*(sigma_Z_res),
                                           ymax = Ztrend + 2*(sigma_Z_res)),
                show.legend = FALSE,
                alpha = 0.2) +  
    geom_ribbon(data = GWLcurves_data, aes(x = -age, ymin = Ztrend - sigma_Z_res,
                                           ymax = Ztrend + sigma_Z_res),
                show.legend = FALSE,
                alpha = 0.3) +  
    geom_ribbon(data = GWLcurves_data, aes(x = -age, ymin = Zxyt - 2*sqrt(var1.var),
                                           ymax = Zxyt + 2*sqrt(var1.var)),
                fill = "darkorange",
                show.legend = FALSE,
                alpha = 0.3) +
    geom_ribbon(data = GWLcurves_data, aes(x = -age, ymin = Zxyt - sqrt(var1.var),
                                           ymax = Zxyt + sqrt(var1.var)),
                fill = "darkorange",
                show.legend = FALSE,
                alpha = 0.5) +
    geom_line(data = GWLcurves_data, aes(x = -age, y = Zxyt),
              linewidth = 1, 
              linetype = "solid",
              # color = "navyblue") +
              color = "gray10") +
    geom_line(data = GWLcurves_data, aes(x = -age, y = Ztrend),
              linewidth = 1, linetype = "longdash", color = "gray40") +
    geom_crossbar(SLIP_data, mapping = aes(x = -t, y = Z,
                                           ymin = Z - Z_2σ_min, ymax = Z + Z_2σ_plus),
                  color = "black", width = SLIP_data$Age_2σ_min - SLIP_data$Age_2σ_plus, fatten = 1,
                  linewidth = 1.5) +
    geom_linerange(ULD_data, mapping = aes(x = -t,
                                           y = Z + Z_2σ_plus,
                                           ymin = Z - Z_2σ_min,
                                           ymax = Z + Z_2σ_plus, 
                                           color = as.factor(type)),
                   linewidth = 1.5) +
    geom_linerange(ULD_data, mapping = aes(x = -t,
                                           y = Z + Z_2σ_plus,
                                           xmin = -Age_2σ_min,
                                           xmax = -Age_2σ_plus, color = as.factor(type)),
                   linewidth = 1.5) +
    scale_color_manual(name = "Indicator type", values = type_colors, labels = c("tidal-ULD", "river gradient ULD", "local ULD")) +
    scale_x_continuous(breaks = c(-12000, -9000, -6000, -3000, 0),  
                       labels = c(12000, 9000, 6000, 3000, 0)) +  
    geom_label(aes(x = -tmin, y = 0, label = run_info),  
               fill = "white", color = "black", fontface = "bold.italic", hjust = 0) + 
    labs(x = "Age [cal. years BP]", 
         y = "Depth [m]") + 
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      legend.position = "none"
    ) +
    ggtitle(paste(name, "output after kriging"))
  
}

############# Store data in Stars format for further analysis #############

star_all <- input_to_stars(total_naomit, FALSE, c(4,5,6,7,11,12,13), 
                           1000, 1000)

# Change all artificial zero values into NA (so they won't be plotted)
star_all[star_all == 0 ]  <- NA

# create stars array with ages as labels, choose filtered below PDEM or all data
time_extract = star_all[,,,]

# Create a label array with all the timestep ages and add "BP" behind it.
age_names <- paste(-t, "BP", sep = " ")

time_extract =  st_set_dimensions(time_extract, "age", values = age_names)
time_extract =  st_set_dimensions(time_extract, "age", values = -t)

timestep1000 = seq(5,length(t),5)
t[timestep1000]


################### Prepare observations to enable facet wrap #############
round_to_200 <- function(x) {
  round(x / 200) * 200
}

round_to_1k <- function(x) {
  round(x / 1000) * 1000
}

obs_location <- obs %>%
  select(x, y, t) %>%
  mutate(age_round = round_to_200(t),
         age = round_to_1k(t))


########################### Visualization #################################
PDEM_map <- ggplot() + 
  geom_sf(data = provincies_shape, color=alpha("black",0.2), fill = "grey", alpha = 1) +
  geom_stars(data = time_extract[2,,,timestep1000], alpha = 1) +
  scale_fill_gradientn(colours = viridis(256, option = "D"),
                       name = "PDEM [m]",
                       breaks = c(-35, -10, 0, 10, 20),
                       labels = c("-35", "-10", "0", "10", "20"),
                       na.value = NA)+
  ggtitle("Extended PDEM") +
  geom_sf(data = provincies_shape, color=alpha("black",0.2), fill = NA, alpha = 0.5) +
  geom_sf(data = Belgium, color=alpha("black",1), fill = "grey", alpha = 1)+
  geom_sf(data = Germany, color=alpha("black",1), fill = "grey", alpha = 1)+
  geom_sf(data = provincies_shape, color=alpha("black",0.2), fill = NA, alpha = 0.5) +
  geom_point(aes(x = obs$x, y = obs$y), 
             size = 2, shape = 16, colour = alpha("black",0.3)) +
  coord_sf(xlim = c(6000, 270000), ylim = c(365000, 630000)) +
  guides(fill = guide_colourbar(title.position = "right"))+
  theme(axis.text.x = element_blank(),
        # theme(axis.text.x = element_text(size = 10, angle = 90),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_blank(),
        # theme(axis.text.x = element_text(size = 10, angle = 90),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 30, face = "bold")) +
  theme(legend.position = "right", legend.key.height = unit(2, "cm"),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30, angle = -90, face = "bold"),
        legend.title.align = 0.5,
        legend.direction = "vertical",
        strip.text.x = element_text(size = 20, color = "black", face = "bold.italic"),
        plot.title = element_text(hjust = 0,size = 20, face = "bold"))


PDEM_map



Zn_map <- plot_interpolation_maps(1, "trend Zn", timestep1000, obs_location)
PDEM_map1 <- plot_interpolation_maps(2, "PDEM [m]", 1, obs_location)
VLM_map <- plot_interpolation_maps(3, "VLM rate [m/kyr]", 1, obs_location)
Z_trend <- plot_interpolation_maps(4, "Trend prediction [m]", timestep1000, obs_location)
Kriging_pred <- plot_interpolation_maps(5, "Kriging prediction [m]", timestep1000, obs_location)
Kriging_var <- plot_interpolation_maps(6, "Kriging variance", timestep1000, obs_location)
Zxyt <- plot_interpolation_maps(7, "Total Zxyt [m]", timestep1000, obs_location)

Zn_map_filename <- save_as_png(output_loc, indicator, correction, fitmodel_name, "Zn_maps")
PDEM_map1_filename <- save_as_png(output_loc, indicator, correction, fitmodel_name, "PDEM_map")
VLM_map_filename <- save_as_png(output_loc, indicator, correction, fitmodel_name, "VLM_rate_map")
Z_trend_filename <- save_as_png(output_loc, indicator, correction, fitmodel_name, "Ztrend_maps")
Kriging_pred_filename <- save_as_png(output_loc, indicator, correction, fitmodel_name, "KrigePred_maps")
Kriging_var_filename <- save_as_png(output_loc,indicator, correction, fitmodel_name, "KrigeVar_maps")
Zxyt_filename <- save_as_png(output_loc, indicator, correction, fitmodel_name,  "Zxyt_maps")


save_map_plots(Zn_map_filename, 2500, 1000, Zn_map)
save_map_plots(PDEM_map1_filename, 1000, 800, PDEM_map1)
save_map_plots(VLM_map_filename, 1000, 800, VLM_map)
save_map_plots(Z_trend_filename, 2500, 1000, Z_trend)
save_map_plots(Kriging_pred_filename, 2500, 1000, Kriging_pred)
save_map_plots(Kriging_var_filename, 2500, 1000, Kriging_var)
save_map_plots(Zxyt_filename, 2500, 1000, Zxyt)


# Standard age-depth plots -------------------------------------------------
Hillegersberg_plot <- ad_plot_single(93600, 441000, "Hillegersberg")
Uitgeest_plot <- ad_plot_single(110000, 504000, "Uitgeest")
BlauweStenen_plot <- ad_plot_single(181000, 520000, "Blauwe Stenen")
Mieden_plot <- ad_plot_single(209000, 585000, "Mieden")
Kappersbult_plot <- ad_plot_single(238000, 572000, "Kappersbult")
Hefswal_plot <- ad_plot_single(245000, 607000, "Hefswal")
Winschoten_plot <- ad_plot_single(269000, 578000, "Winschoten")
Koegras_plot <- ad_plot_single(112000, 548000, "Koegras")
Leerdam_plot <- ad_plot_single(135000, 436000, "Leerdam")
OudeStoof_plot <- ad_plot_single(56000, 374000, "Oude Stoof")

grid.arrange(Hefswal_plot, Mieden_plot,  
             Winschoten_plot, Kappersbult_plot,
             nrow = 2)

grid.arrange(Koegras_plot, Uitgeest_plot, BlauweStenen_plot, 
             Hillegersberg_plot, Leerdam_plot, OudeStoof_plot,
             nrow = 2)


save_ad_plots <-  function(filename, width, height, plot_name) {
  png(file = filename, width = width, height = height)
  print({
    plot_name
  })
  
  dev.off()
  
}


singleplots_1_filename <- save_as_png(output_loc, indicator, correction, fitmodel_name, "singleADplots_1")
singleplots_2_filename <- save_as_png(output_loc, indicator, correction, fitmodel_name, "singleADplots_2")

save_ad_plots(singleplots_1_filename, 1500, 1000, 
              grid.arrange(Hefswal_plot, Mieden_plot,  
                           Winschoten_plot, Kappersbult_plot,
                           nrow = 2))
save_ad_plots(singleplots_2_filename, 1500, 1000,  
              grid.arrange(Koegras_plot, Uitgeest_plot, BlauweStenen_plot, 
                           Hillegersberg_plot, Leerdam_plot, OudeStoof_plot,
                           nrow = 2))
