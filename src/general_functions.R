# Kim de Wit
# The code  presented in this script is part of the 
# “Living on Soft Soils: Subsidence and Society” project 
# (grant no. NWA.1160.18.259), WP1.3

###########################################################################
# General functions for paleo-water level interpolation
###########################################################################
# Script containing general functions used in the paleo-water level 
# interpolation. Required for running other interpolation scripts

################################ Libraries ################################

library(readxl)
library(snakecase)
library(gstat)

############################### Functions #################################

adjust_headers <- function (data) {
  
  headers_snake <- to_snake_case(colnames(data))
  headers_snake[8] <- c("region_name")
  headers_snake[16:18] <- c("age_cal_a_bp",
                            "age_2_σ_uncertainty_cal_a_plus",
                            "age_2_σ_uncertainty_cal_a_min")
  headers_snake[20:21] <- c("ox_cal_modelled_age_2_σ_uncertainty_cal_a_plus",
                            "ox_cal_modelled_age_2_σ_uncertainty_cal_a_min")
  headers_snake[50:51] <- c("sample_elevation_uncertainty_m_plus",
                            "sample_elevation_uncertainty_m_min")
  headers_snake[81:82] <- c("z_gw_uncertainty_m_plus",
                            "z_gw_uncertainty_m_min")
  headers_snake[84:85]  <- c("rsl_2_σ_uncertainty_m_plus",
                             "rsl_2_σ_uncertainty_m_min")
  headers_snake[87:88] <- c("corrected_gwl_2_σ_uncertainty_m_plus",
                            "corrected_gwl_2_σ_uncertainty_m_min")
  headers_snake[90:91] <- c("corrected_rsl_2_σ_uncertainty_m_plus",
                            "corrected_rsl_2_σ_uncertainty_m_min")
  headers_snake[94:95] <- c("withBG_tect_corrected_rsl_2_σ_uncertainty_m_plus",
                            "withBG_tect_corrected_rsl_2_σ_uncertainty_m_min")
  
  colnames(data) <- headers_snake
  colnames(data)[1] == "unique_sample_id"
  colnames(data)[1] <- "uniqueID"
  
  data

}


create_raster <- function(filelocation, crs){
  #' Create rasters
  #' 
  #' @description Creates rasters with desired coordinate projection from 
  #' input file
  #' @return Raster with assigned coordinate projection
  #' @details The input should be a raster file (e.g. .tiff)
  raster <- read_stars(filelocation)
  raster_with_crs <- st_set_crs(raster, crs)
  return(raster_with_crs)
}


input_to_stars <- function (input_data, above_PDEM, attributes, dx, dy) {

  if (above_PDEM == TRUE) {
      data <- na.omit(input_data)
      data[data$Zxyt < data$PDEM, 13] = NA
      data[data$Ztrend < data$PDEM, 7] = NA
  } else {
    data <- na.omit(input_data)
  }
  

  for (j in attributes){

    # Set input data frame for stars
    inputdf =  data.frame(att = data[,j], age = data$age, x = data$x, 
                          y = data$y)
    
    # Create one list with all lists/stars of each time step 
    C_test = lapply(-t, function(x) {
      inputdf[inputdf$age == x, ]
    })
    
    # Create list with rasters from each time step list
    total_raster  = lapply(seq(1,length(t),1), function(x) {
      listtmp <- st_as_sf(C_test[[x]], coords = c("x","y","age"), remove = TRUE, crs = RDnew)
      st_rasterize(listtmp, dx = dx, dy = dy)
    })
    
    # Asign the time steps as name to each individual raster
    names(total_raster) <- -t
    
    # Turn the rasters into one stars object and again asign correct names
    star = do.call("c", total_raster)
    names(star) <- t
    
    # Merge the timesteps into a new dimension "age" of the stars
    star = merge(star, name = "age")
    
    # Rename the attribute to the original attribute name
    names(star) = colnames(data[j])
    
    if (j == 4) {
      # Add initial attribute
      output_star = c(star)  
    } else {
      # Stack the stars into one star with multiple attributes
      output_star = c(output_star, star)
    } 
  }
  
  # Change all artificial zero values into NA (so they won't be plotted)
  output_star[output_star == 0 ]  <- NA
  
  output_star
}


filter_output_slope <- function(input_data, slope_lim, age, outputname){
  input_split <- split(input_data)
  input_split_df <- as.data.frame(input_split)
  
  slope_data <- starsExtra::slope(input_split)
  slope_data$age <- age
  slope_data$slope <- units::drop_units(slope_data$slope)
  slope_df <- as.data.frame((slope_data))
  
  slope_df$data <- input_split_df[,3]
  slope_df <- slope_df %>%
    mutate(
      # Set data to NA if slope exceeds slope_lim
      data = if_else((!is.na(slope) & slope >= slope_lim), 
                     NA_real_, data)
    )
  
  colnames(slope_df)[5] <- outputname
  
  slope_df
}


filter_Z_data <- function(input_data, min_rate, timesteps){
  
  # Apply max slope and minimal rise rate filter to output data
  for (timestep in 1:length(timesteps)) { 
    
    diff_list <- vector("list", length = length(timesteps)-1)
    
    for (i in 1:(length(timesteps)-1)) {
      # Subtract the current time step from the previous one and divide by the difference in time to calculate the rate
      # multiple by 1000 to convert rate to [mm/yr]
      diff_list[[i]] <- (input_data[7,,,i+1, drop = TRUE] - input_data[7,,,i, drop = TRUE])/(timesteps[1+i]-timesteps[i])*1000
    }
    
    rates_of_change <- st_as_stars(do.call(c, diff_list))
    
    if (timestep == 1) {
      mask <- rates_of_change[1] > min_rate[timestep]
      
      rate_filtered_Zxyt <- input_data[7,,,timestep]
      rate_filtered_Zxyt[!mask] <- NA
      
      rate_filtered_Ztrend <- input_data[4,,,timestep]
      rate_filtered_Ztrend[!mask] <- NA
      
    } else {
      mask <- rates_of_change[timestep-1] > min_rate[timestep]
      
      rate_filtered_Zxyt <- input_data[7,,,timestep]
      rate_filtered_Zxyt[!mask] <- NA
      
      rate_filtered_Ztrend <- input_data[4,,,timestep]
      rate_filtered_Ztrend[!mask] <- NA
      
    }
    
    filtered_Zxyt <- filter_output_slope(rate_filtered_Zxyt, 0.02, -timesteps[timestep], "Zxyt")
    filtered_Ztrend <- filter_output_slope(rate_filtered_Ztrend, 0.02, -timesteps[timestep], "Ztrend")
    
    filtered_data_tstep <- filtered_Zxyt
    filtered_data_tstep$Ztrend <- filtered_Ztrend$Ztrend
    
    # Add kriging variance
    krige_var_star <- input_data[6,,,]
    # krige_var_star <- interpolation_output[6,,,]
    krige_var_df <- as.data.frame(krige_var_star) 
    krige_var_df <- krige_var_df %>%
      filter(age %in% -timesteps[timestep])
    
    filtered_data_tstep$var <- krige_var_df$var1.var
    
    # Add VLM-rates
    VLM_star <- input_data[3,,,]
    # krige_var_star <- interpolation_output[6,,,]
    VLM_df <- as.data.frame(VLM_star) 
    VLM_df <- VLM_df %>%
      filter(age %in% -timesteps[timestep])
    
    filtered_data_tstep$VLM_rate <- VLM_df$VLM_rate
    
  
    
    
    if (timestep == 1) {
      filtered_data = filtered_data_tstep
    } else {
      filtered_data <- rbind(filtered_data, filtered_data_tstep)
    }
  }
  
  filtered_data
  
}

