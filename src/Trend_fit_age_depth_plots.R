# Kim de Wit
# The code  presented in this script is part of the 
# “Living on Soft Soils: Subsidence and Society” project 
# (grant no. NWA.1160.18.259), WP1.3

###########################################################################
# Trend fit visualization script
###########################################################################
# Script for visualizing age-dept plots at specific location, useful for 
# analyzing the output of the fitted trend.
# This script should be ran directly after "Trend-fit-script.R", to ensure 
# all required input data is loaded.

################################ Libraries ################################
library(cowplot)
library(ggrepel)


### Plot functions ###############
ad_plot_single <-  function(x, y, name) {
  plot_loc <- data.frame(x = x, y = y, location = name)
  
  # Extract GWL1000, GWL10800, Dxy and q at observation locations
  plot_loc_points <- st_as_sf(plot_loc, coords = c("x","y"), crs = RDnew)
  
  plot_loc_extract = data.frame(ZA0 = st_extract(ZA0_raster, plot_loc_points),
                                ZA1 = st_extract(ZA1_raster, plot_loc_points),
                                Dxy = st_extract(Dxy_raster, plot_loc_points),
                                q = st_extract(q_raster, plot_loc_points))
  
  # Transform data
  plot_loc$x_n = 1 -  (mean(grid2D$x) -  plot_loc$x)/(xmax-xmin)
  plot_loc$y_n = 1 -  (mean(grid2D$y) -  plot_loc$y)/(ymax-ymin)
  
  # Calculate parameters
  plot_loc$a <- calc_a(params, plot_loc$x_n, plot_loc$y_n)
  plot_loc$c <- calc_c(params, plot_loc$x_n, plot_loc$y_n)
  plot_loc$b <- params[index_b]
  
  # At extracted DEM data to observations
  plot_loc$ZA0 = plot_loc_extract[[1]]
  plot_loc$ZA1 = plot_loc_extract[[3]]
  plot_loc$Dxy = plot_loc_extract[[5]]
  plot_loc$q = plot_loc_extract[[7]]
  
  locs_plotdata <- data.frame()
  for(location in 1:length(plot_loc$location)) {
    plotdata <- data.frame(location = plot_loc$location[location],
                           x = plot_loc$x[location],
                           y = plot_loc$y[location],
                           age = t,
                           Zn = sigmoid_function(plot_loc$a[location], 
                                                 plot_loc$b[location], 
                                                 plot_loc$c[location], plot_loc$q[location], p, 
                                                 plot_loc$x_n[location], plot_loc$y_n[location])
    )
    
    plotdata$Z <- plot_loc$ZA0[location] + plotdata$Zn * plot_loc$Dxy[location]
    locs_plotdata <- rbind(locs_plotdata, plotdata)
  }
  
  trend <- locs_plotdata
  
  obs_loc <- obs %>%
    filter(grepl(name, Name))
  
  SLIP_data <- obs_loc %>%
    filter(type %in% "0")
  
  ULD_data <- obs_loc %>%
    filter(type %in% c("1", "2", "3"))
  
  run_info = paste(region_name, indicator, "data", correction, "\n trend fit:", 
                   fitmodel_name, "\n Sigma:", round(sigma_Z, 2), "m", "\n n:", n, sep = " ")
  
  type_colors <- c("1" = "#0072B2", "2" = "lightseagreen", "3" = "#D55E00")
  
  ggplot() +
    geom_ribbon(data = trend, aes(x = age, ymin = Z - sigma_Z,
                                  ymax = Z + sigma_Z),
                show.legend = FALSE,
                alpha = 0.2) +  # Ribbon for SD
    geom_line(data = trend, aes(x = age, y = Z), 
              linewidth = 2) +
    geom_crossbar(SLIP_data, mapping = aes(x = -t, y = Z, 
                                           ymin = Z - Z_2σ_min, ymax = Z + Z_2σ_plus), 
                  color = "red", width = SLIP_data$Age_2σ_min - SLIP_data$Age_2σ_plus, fatten = 1) +
    geom_crossbar(SLIP_data, mapping = aes(x = -t, y = Z, 
                                           ymin = Z - Z_2σ_min, ymax = Z + Z_2σ_plus), 
                  color = "red", width = SLIP_data$Age_2σ_min - SLIP_data$Age_2σ_plus, fatten = 1) +
    
    geom_linerange(ULD_data, mapping = aes(x = -t, 
                                           y = Z + Z_2σ_plus, 
                                           ymin = Z - Z_2σ_min, 
                                           ymax = Z + Z_2σ_plus, color = as.factor(type)),
                   linewidth = 1) + 
    geom_linerange(ULD_data, mapping = aes(x = -t, 
                                           y = Z + Z_2σ_plus, 
                                           xmin = -Age_2σ_min, 
                                           xmax = -Age_2σ_plus, color = as.factor(type)),
                   linewidth = 1) + 
    scale_color_manual(name = "Indicator type", values = type_colors, labels = c("tidal-ULD", "river gradient ULD", "local ULD")) +
    scale_x_continuous(breaks = c(-12000, -9000, -6000, -3000, 0),  
                       labels = c(12000, 9000, 6000, 3000, 0)) +  
    geom_label(aes(x = -tmin, y = 0, label = run_info), 
               fill = "white", color = "black", fontface = "bold.italic", hjust = 0) + 
    labs(x = "Age [cal. years BP]", 
         y = "Depth [m]") + 
    theme_bw()+
    theme(
      legend.title = element_text(colour = "black", size = 15, face ="bold"),
      legend.text = element_text(colour = "black", size = 14),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      legend.position = c(.98, .02),
      legend.justification = c("right", "bottom"),
      legend.box.just = "right",
      legend.margin = margin(2, 6, 6, 6),
      legend.background = element_rect(
        fill = alpha("white", 1),
        color = "black",
        linewidth = 0.5
      )
    ) +
    ggtitle(paste(name, "trend before kriging"))
  
}

ad_plot_multiple <-  function(x, y, names, groupname, ribbon) {
  plot_loc <- data.frame(x = x, y = y, location = names)
  
  # Extract GWL1000, GWL10800, Dxy and q at observation locations
  plot_loc_points <- st_as_sf(plot_loc, coords = c("x","y"), crs = RDnew)
  
  plot_loc_extract = data.frame(ZA0 = st_extract(ZA0_raster, plot_loc_points),
                                ZA1 = st_extract(ZA1_raster, plot_loc_points),
                                Dxy = st_extract(Dxy_raster, plot_loc_points),
                                q = st_extract(q_raster, plot_loc_points))
  
  # Transform data
  plot_loc$x_n = 1 -  (mean(grid2D$x) -  plot_loc$x)/(xmax-xmin)
  plot_loc$y_n = 1 -  (mean(grid2D$y) -  plot_loc$y)/(ymax-ymin)
  
  # Calculate parameters
  plot_loc$a <- calc_a(params, plot_loc$x_n, plot_loc$y_n)
  plot_loc$c <- calc_c(params, plot_loc$x_n, plot_loc$y_n)
  plot_loc$b <- params[index_b]
  
  # Obtain ZA0 and Dxy from input raster
  plot_loc$ZA0 = plot_loc_extract[[1]]
  plot_loc$ZA1 = plot_loc_extract[[3]]
  
  # TODO check what is meant here. Is this correct?
  # Using LDEM as ZA0
  plot_loc$Dxy = plot_loc_extract[[5]]
  plot_loc$q = plot_loc_extract[[7]]
  
  locs_plotdata <- data.frame()
  for(location in 1:length(plot_loc$location)) {
    plotdata <- data.frame(location = plot_loc$location[location],
                           x = plot_loc$x[location],
                           y = plot_loc$y[location],
                           age = t,
                           Zn = sigmoid_function(plot_loc$a[location], 
                                                 plot_loc$b[location], 
                                                 plot_loc$c[location], plot_loc$q[location], p, 
                                                 plot_loc$x_n[location], plot_loc$y_n[location])
    )
    
    plotdata$Z <- plot_loc$ZA0[location] + plotdata$Zn * plot_loc$Dxy[location]
    locs_plotdata <- rbind(locs_plotdata, plotdata)
  }
  
  run_info = paste(region_name, indicator, "data", correction, "\n trend fit:", 
                   fitmodel_name, "\n Sigma:", round(sigma_Z, 2), "m", "\n n:", n, sep = " ")
  
  if(ribbon == TRUE) {
    ggplot() +
      geom_ribbon(data = locs_plotdata, aes(x = age, ymin = Z - sigma_Z,
                                            ymax = Z + sigma_Z, fill = location),
                  show.legend = FALSE,
                  alpha = 0.2) +  # Ribbon for SD
      geom_line(data = locs_plotdata, aes(x = age, y = Z, color = location, 
                                          linetype = as.factor(y)), 
                linewidth = 2) +
      scale_fill_viridis_d(option = "D") +
      scale_color_viridis_d(option = "D") +
      scale_x_continuous(breaks = c(-9000, -6000, -3000, 0),  # Customize the breaks (where the ticks appear)
                         labels = c(9000, 6000, 3000, 0)) +  # Customize the labels
      geom_label(aes(x = -tmin, y = 0, label = run_info),  # Add a label at a specific point
                 fill = "white", color = "black", fontface = "bold.italic", hjust = 0) +  # Text box properties
      labs(x = "Age [cal. years BP]", 
           y = "Depth [m]",
           color = "Location",
           linetype = "y-coordinate") + 
      ggtitle(paste(groupname, "trend before kriging")) + 
      theme_bw()+
      theme(
        legend.title = element_text(colour="black", size = 15, face="bold"),
        legend.text = element_text(colour="black", size = 14),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.title = element_text(size = 14, face="bold"),
        axis.text = element_text(size = 12, face="bold"),
        legend.position = c(.98, .02),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(2, 6, 6, 6),
        legend.background = element_rect(
          fill = alpha("white", 1),
          color = "black",
          linewidth = 0.5
        )
      ) +  
      guides(linetype = guide_legend(override.aes = list(size = 0.5),   # Make legend lines smaller
                                     keywidth = 3.5, keyheight = 1.2))    # Adjust legend line width
  } else{
    ggplot() +
      geom_line(data = locs_plotdata, aes(x = age, y = Z, color = location, 
                                          linetype = as.factor(y)), 
                linewidth = 2) +
      scale_color_viridis_d(option = "D") +
      scale_x_continuous(breaks = c(-9000, -6000, -3000, 0),  # Customize the breaks (where the ticks appear)
                         labels = c(9000, 6000, 3000, 0)) +  # Customize the labels
      # geom_label(aes(x = -tmin, y = 0, label = run_info),  # Add a label at a specific point
      #            fill = "white", color = "black", fontface = "bold.italic", hjust = 0) +  # Text box properties
      labs(x = "Age [cal. years BP]", 
           y = "Depth [m]",
           color = "Location",
           linetype = "y-coordinate")+ 
      ggtitle(paste(groupname, "trend before kriging")) + 
      theme_bw()+
      theme(
        legend.title = element_text(colour="black", size = 15, face="bold"),
        legend.text = element_text(colour="black", size = 14),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.title = element_text(size = 14, face="bold"),
        axis.text = element_text(size = 12, face="bold"),
        legend.position = c(.98, .02),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(2, 6, 6, 6),
        legend.background = element_rect(
          fill = alpha("white", 1),
          color = "black",
          linewidth = 0.5
        )
      ) +  
      guides(linetype = guide_legend(override.aes = list(size = 0.5),   # Make legend lines smaller
                                     keywidth = 3.5, keyheight = 1.2))    # Adjust legend line width
    
    
    
  }
  
  
}

map_loc_multiple <-  function(x, y, name) {
  plot_loc <- data.frame(x = x, y = y, location = name)
  
  ggplot() +  
    geom_sf(data = provincies_shape, color = alpha("black",0.2), fill = "grey", alpha = 1) +
    geom_sf(data = studyarea_shape, color = alpha("black", 0.8), alpha = 0, 
            linetype = "dashed") +  
    geom_point(data = plot_loc, aes(x, y, color = location), size = 3) +
    scale_color_viridis_d(option = "D") +
    geom_text_repel(data = plot_loc, aes(x, y, label = location), 
                    box.padding = 0.5, # Adjust space between labels and points
                    max.overlaps = Inf  # Allows all labels to be displayed
    ) +
  theme(axis.title = element_blank(),  # Removes axis titles
        axis.text = element_blank(),   # Removes axis labels
        axis.ticks = element_blank(),   # Removes axis ticks
        legend.position = "none",      # Removes the legend
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
    coord_sf(xlim = c(6000, 270000), ylim = c(355000, 630000))

  
  
}

save_ad_plots <-  function(filename, width, height, plot_name) {
  png(file = filename, width = width, height = height)
  print({
    plot_name
  })
  
  dev.off()
  
}

# Function for saving as svg
# file location
save_as_svg <- function(fileloc_base, region_name, indicator, correction,
                        fitmodel_name, variable){ 
  paste(fileloc_base, region_name, indicator, correction, fitmodel_name,
        paste0(fitmodel_name, region_name, "_", variable, ".svg"),sep = "/")
}

save_ad_plots_svg <-  function(filename, width, height, plot_name) {
  
  svg(file = filename,
      width = width, height = height) 
  
  plot(plot_name)
  
  # Close the device
  dev.off()
  
}

map_Zresiduals_timesliced <- function(data, lower_limit, upper_limit) {
  timeslice_data <- trendData %>%
    filter(t >= lower_limit & t <= upper_limit)
  
  ggplot() +  
    geom_sf(data = provincies_shape, color = alpha("black",0.6), fill = "grey88", alpha = 1) +
    geom_point(aes(x = timeslice_data$x, y = timeslice_data$y, 
                   fill = timeslice_data$resZ, 
                   size = abs(timeslice_data$resZ)), 
               shape = 21, colour = "black", stroke = 0.5) +  
    geom_label(aes(x = 20000, y = 600000, 
                   label = paste("Map of Z residuals \n", "n = ", length(timeslice_data$x))),  
               fill = "white", color = "black",
               fontface = "bold.italic", hjust = 0) + 
    scale_fill_gradientn(colors = colors_res,
                         name = "Residuals [m]",
                         limits = c(min(trendData$resZ), max(trendData$resZ))) +
    scale_size(breaks = seq(0, 4, 0.5)) +
    guides(size = "none") +
    theme(
      plot.title = element_text(hjust = 0.5),  
      axis.title = element_blank(),    
      axis.text = element_blank(),     
      axis.ticks = element_blank(),     
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", color = "black")  
    ) + 
    guides(fill = guide_colorbar(
      title.position = "top", title.hjust = 0,  
      barwidth = unit(0.5, "cm"), barheight = unit(3, "cm")  
    )) +
    ggtitle(paste(lower_limit, "-", upper_limit, "cal. years BP")) +
    coord_sf(xlim = c(15000, 270000), ylim = c(365000, 610000))
  
}
map_Znresiduals_timesliced <- function(data, lower_limit, upper_limit) {
  timeslice_data <- trendData %>%
    filter(t >= lower_limit & t <= upper_limit)
  
  ggplot() +  
    geom_sf(data = provincies_shape, color = alpha("black",0.6), fill = "grey88", alpha = 1) +
    geom_point(aes(x = timeslice_data$x, y = timeslice_data$y, 
                   fill = timeslice_data$resZn, 
                   size = abs(timeslice_data$resZn)), 
               shape = 21, colour = "black", stroke = 0.5) +  
    geom_label(aes(x = 20000, y = 600000, 
                   label = paste("Map of Zn residuals \n", "n = ", length(timeslice_data$x))), 
               fill = "white", color = "black",
               fontface = "bold.italic", hjust = 0) +  
    scale_fill_gradientn(colors = colors_res,
                         name = "Residuals [m]",
                         limits = c(min(trendData$resZn), max(trendData$resZn))) +
    scale_size(breaks = seq(0, 4, 0.5)) +
    guides(size = "none") +
    theme(
      plot.title = element_text(hjust = 0.5),  
      axis.title = element_blank(),    
      axis.text = element_blank(),     
      axis.ticks = element_blank(),     
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", color = "black")  
    ) + 
    guides(fill = guide_colorbar(
      title.position = "top", title.hjust = 0, 
      barwidth = unit(0.5, "cm"), barheight = unit(3, "cm")  
    )) +
    ggtitle(paste(lower_limit, "-", upper_limit, "cal. years BP")) +
    coord_sf(xlim = c(15000, 270000), ylim = c(365000, 610000))
  
}

### Create multiple locations plots ###############
# Thick Holocene wedge group
ad_thickHolocene <- ad_plot_multiple(c(245000, 177000, 106000, 129000, 93600, 57000), 
                 c(607000, 607000, 524000, 484000, 441000, 408000), 
                 c("Hefswal", "Molengat", "Schoorl", 
                   "Diemen", "Hillegersberg", "Bouwlust"), 
                 "Thick Holocene wedge", FALSE)

ad_thickHolocene_ribbon <- ad_plot_multiple(c(245000, 177000, 106000, 129000, 93600, 57000), 
                                     c(607000, 607000, 524000, 484000, 441000, 408000), 
                                     c("Hefswal", "Molengat", "Schoorl", 
                                       "Diemen", "Hillegersberg", "Bouwlust"), 
                                     "Thick Holocene wedge", TRUE)

map_thickHolocene <- map_loc_multiple(c(245000, 177000, 106000, 129000, 93600, 57000), 
                                      c(607000, 607000, 524000, 484000, 441000, 408000), 
                                      c("Hefswal", "Molengat", "Schoorl", 
                                        "Diemen", "Hillegersberg", "Bouwlust"))

# Inland points group
ad_inland <- ad_plot_multiple(c(261000, 165000, 172000, 112000, 135000, 57000),
                 c(585000, 574000, 510000, 548000, 436000, 374000),
                 c("Nieuwolda", "Winsum", "Swifterbant", 
                   "Koegras", "Leerdam", "Oude Stoof"), 
                 "Inland", FALSE) 

ad_inland_ribbon <- ad_plot_multiple(c(261000, 165000, 172000, 112000, 135000, 57000),
                              c(585000, 574000, 510000, 548000, 436000, 374000),
                              c("Nieuwolda", "Winsum", "Swifterbant", 
                                "Koegras", "Leerdam", "Oude Stoof"), 
                              "Inland", TRUE) 

map_inland <- map_loc_multiple(c(261000, 165000, 172000, 112000, 135000, 57000),
                                 c(585000, 574000, 510000, 548000, 436000, 374000),
                                 c("Nieuwolda", "Winsum", "Swifterbant", 
                                   "Koegras", "Leerdam", "Oude Stoof"))

# Along one y-coordinate
ad_westeast <- ad_plot_multiple(c(57000, 93600, 126000, 159000, 188000, 205000),
                 c(430000, 430000, 430000, 430000, 430000, 430000),
                 c("Coast", "Rotterdam", "Gorinchem", 
                   "Tiel", "Nijmegen", "Lobith"),
                 "West-east", FALSE)

ad_westeast_ribbon <- ad_plot_multiple(c(57000, 93600, 126000, 159000, 188000, 205000),
                                c(430000, 430000, 430000, 430000, 430000, 430000),
                                c("Coast", "Rotterdam", "Gorinchem", 
                                  "Tiel", "Nijmegen", "Lobith"),
                                "West-east", TRUE)

map_westeast <- map_loc_multiple(c(57000, 93600, 126000, 159000, 188000, 205000),
                 c(430000, 430000, 430000, 430000, 430000, 430000),
                 c("Coast", "Rotterdam", "Gorinchem", 
                   "Tiel", "Nijmegen", "Lobith"))


# Combine the two plots
westeast_plot <- ggdraw() +
  draw_plot(ad_westeast) +  
  draw_plot(map_westeast, x = 0.55, y = -0.001, width = 0.25, height = 0.5)
# print(westeast_plot)

thickHolocene_plot <- ggdraw() +
  draw_plot(ad_thickHolocene) +  
  draw_plot(map_thickHolocene, x = 0.55, y = -0.001, width = 0.3, height = 0.5)
# print(thickHolocene_plot)

inland_plot <- ggdraw() +
  draw_plot(ad_inland) +  
  draw_plot(map_inland, x = 0.55, y = -0.001, width = 0.25, height = 0.5)
# print(inland_plot)

# With ribbons
westeast_plot_ribbon <- ggdraw() +
  draw_plot(ad_westeast_ribbon) +  
  draw_plot(map_westeast, x = 0.55, y = -0.001, width = 0.25, height = 0.5)

thickHolocene_plot_ribbon <- ggdraw() +
  draw_plot(ad_thickHolocene_ribbon) +  
  draw_plot(map_thickHolocene, x = 0.55, y = -0.001, width = 0.3, height = 0.5)

inland_plot_ribbon <- ggdraw() +
  draw_plot(ad_inland_ribbon) +  
  draw_plot(map_inland, x = 0.55, y = -0.001, width = 0.25, height = 0.5)


# Combine the two plots
westeast_plot <- ad_westeast

thickHolocene_plot <- ad_thickHolocene

inland_plot <- ad_inland

# With ribbons
westeast_plot_ribbon <- ad_westeast_ribbon

thickHolocene_plot_ribbon <- ad_thickHolocene_ribbon

inland_plot_ribbon <- ad_inland_ribbon


# Optional possibility for saving plots as png-files
# # Specify output location
# ad_inland_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "ad_Inland")
# ad_thick_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "ad_ThickHolocene")
# ad_westeast_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "ad_westeast")
# ad_inland_filename_ribbon <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "ad_Inland_ribbon")
# ad_thick_filename_ribbon <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "ad_ThickHolocene_ribbon")
# ad_westeast_filename_ribbon <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "ad_westeast_ribbon")
# 
# # Save plots
# save_ad_plots(ad_inland_filename, 1000, 900, inland_plot)
# save_ad_plots(ad_thick_filename, 1000, 900, thickHolocene_plot) 
# save_ad_plots(ad_westeast_filename, 1000, 900, westeast_plot) 
# save_ad_plots(ad_inland_filename_ribbon, 1000, 900, inland_plot_ribbon)
# save_ad_plots(ad_thick_filename_ribbon, 1000, 900, thickHolocene_plot_ribbon) 
# save_ad_plots(ad_westeast_filename_ribbon, 1000, 900, westeast_plot_ribbon) 


# svg-plots
# Specify output location
ad_inland_filename <- save_as_svg(output_loc, region_name, indicator, correction, fitmodel_name, "ad_Inland")
ad_thick_filename <- save_as_svg(output_loc, region_name, indicator, correction, fitmodel_name, "ad_ThickHolocene")
ad_westeast_filename <- save_as_svg(output_loc, region_name, indicator, correction, fitmodel_name, "ad_westeast")
ad_inland_filename_ribbon <- save_as_svg(output_loc, region_name, indicator, correction, fitmodel_name, "ad_Inland_ribbon")
ad_thick_filename_ribbon <- save_as_svg(output_loc, region_name, indicator, correction, fitmodel_name, "ad_ThickHolocene_ribbon")
ad_westeast_filename_ribbon <- save_as_svg(output_loc, region_name, indicator, correction, fitmodel_name, "ad_westeast_ribbon")

# Save plots
save_ad_plots_svg(ad_inland_filename, 10, 9, inland_plot)
save_ad_plots_svg(ad_thick_filename, 10, 9, thickHolocene_plot) 
save_ad_plots_svg(ad_westeast_filename, 10, 9, westeast_plot) 
save_ad_plots_svg(ad_inland_filename_ribbon, 10, 9, inland_plot_ribbon)
save_ad_plots_svg(ad_thick_filename_ribbon, 10, 9, thickHolocene_plot_ribbon) 
save_ad_plots_svg(ad_westeast_filename_ribbon, 100, 9, westeast_plot_ribbon) 


### Create Plots of single locations #######################################
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



map_single_plots <- map_loc_multiple(c(93600, 110000, 181000, 209000, 238000, 
                                       245000, 269000, 112000, 135000, 56000),
                                     c(441000, 504000, 520000, 585000, 572000, 
                                       607000, 578000, 548000, 436000, 374000),
                                     c("Hillegersberg", "Uitgeest", "Blauwe Stenen", 
                                       "Mieden", "Kappersbult", "Hefswal", "Winschoten", 
                                       "Koegras", "Leerdam", "Oude Stoof"))

ad_Hillegersberg_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "ad_Hillegersberg")
save_ad_plots(ad_Hillegersberg_filename, 1000, 900, Hillegersberg_plot)

singleplots_1_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "singleADplots_1")
singleplots_2_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "singleADplots_2")
singleplots_map_filename <- save_as_png(output_loc, region_name, indicator, correction, fitmodel_name, "SingleSitesMap")

save_ad_plots(singleplots_1_filename, 1500, 1000, 
              grid.arrange(Hefswal_plot, Mieden_plot,  
                           Winschoten_plot, Kappersbult_plot,
                           nrow = 2))
save_ad_plots(singleplots_2_filename, 1500, 1000,  
              grid.arrange(Koegras_plot, Uitgeest_plot, BlauweStenen_plot, 
                           Hillegersberg_plot, Leerdam_plot, OudeStoof_plot,
                           nrow = 2))
save_ad_plots(singleplots_map_filename, 850, 900, map_single_plots)


### time sliced residuals maps #############################################
Zres_maps_filename <- save_as_png(output_loc, region_name, indicator, 
                                  correction, fitmodel_name, "Zres_timeMaps")
Znres_maps_filename <- save_as_png(output_loc, region_name, indicator, 
                                  correction, fitmodel_name, "Znres_timeMaps")

save_ad_plots(Zres_maps_filename, 2500, 1000, 
              grid.arrange(map_Zresiduals_timesliced(trendData, 1000, 2000),
                           map_Zresiduals_timesliced(trendData, 2000, 3000),
                           map_Zresiduals_timesliced(trendData, 3000, 4000),
                           map_Zresiduals_timesliced(trendData, 4000, 5000),
                           map_Zresiduals_timesliced(trendData, 5000, 6000),
                           map_Zresiduals_timesliced(trendData, 6000, 7000),
                           map_Zresiduals_timesliced(trendData, 7000, 8000),
                           map_Zresiduals_timesliced(trendData, 8000, 9000),
                           map_Zresiduals_timesliced(trendData, 9000, 10000),
                           map_Zresiduals_timesliced(trendData, 10000, 11000), 
                           nrow = 2))

save_ad_plots(Znres_maps_filename, 2500, 1000, 
              grid.arrange(map_Znresiduals_timesliced(trendData, 1000, 2000),
                           map_Znresiduals_timesliced(trendData, 2000, 3000),
                           map_Znresiduals_timesliced(trendData, 3000, 4000),
                           map_Znresiduals_timesliced(trendData, 4000, 5000),
                           map_Znresiduals_timesliced(trendData, 5000, 6000),
                           map_Znresiduals_timesliced(trendData, 6000, 7000),
                           map_Znresiduals_timesliced(trendData, 7000, 8000),
                           map_Znresiduals_timesliced(trendData, 8000, 9000),
                           map_Znresiduals_timesliced(trendData, 9000, 10000),
                           map_Znresiduals_timesliced(trendData, 10000, 11000), 
                           nrow = 2))

