#--------------------------------------------------------------------------------------
# Hillslope wizzard developed by Ralf Loritz (2015)
# Version 0.1
# Delination of 2D catenas from hillslopes. The hillslope, dist2river, elev2river and accumulation maps
# where extracte using Whitebox GIS in combination Grass GIS. The idea behind this code is based on a approach from 
# Francke et al. 2006 (Automated catena-based discretization of landscapes for the derivation of hydrological modelling units)
# and Cochrane and Flanagan (2003). 
# The idea behind this approach is to select every unique distance and there corrosponding elevation from a river inside a hillslope.
# If more than one raster cell have the same distance the elevation is calculated by a weighted mean using the flow accumulation.
# The hillslope width can be estimated by the frequency of the appearance of a single dinstance. 

#--------------------------------------------------------------------------------------
#Load additional packages (used package version 2.2-31)
library('raster')
library('rgdal') # this could be changed to sf in case of dependency clashes with rgdal
library('dplyr')
library('tidyr')
library('ggplot2')
#--------------------------------------------------------------------------------------
# Own packages
source("plot.R")

#--------------------------------------------------------------------------------------
# Start hillslope method
#--------------------------------------------------------------------------------------
hillslope_tool <- function(hillslope_nr, li_spatial, plot_2d_catena=FALSE, plot_hillslope_width=FALSE,
                           no_rf=0.2, min_dist=3, min_area=10000, freedom=10)
{
  hillslopes <- li_spatial$hillslopes
  dist2river <- li_spatial$dist2river
  elev2river <- li_spatial$elev2river
  aspect <- li_spatial$aspect
  accum <- li_spatial$accum
  stream_id <- li_spatial$stream_id
  hillslope_as_pts <- li_spatial$hillslope_table
  landuse <- li_spatial$landuse
  
  # safty check if hillslope_nr is set to zero
  if(hillslope_nr == 0) stop("---Do not use hillslope_nr = 0---")
  
  if (hillslope_nr == -1) {
    hill <- hillslope_as_pts  # Use all hillslopes when hillslope_nr == -1
  } else {
    hill <- hillslope_as_pts[hillslope_as_pts[, 3] == hillslope_nr, ]
  }
  
  #extract number of cells and area of hillslope
  number_of_cells <- length(hill[,1])
  area <- number_of_cells * res(hillslopes)[1]^2 # number of pixel in a hillslope * resolution of raster
  
  #--------------------------------------------------------------------------------------
  # select hillslope specific elevations, distances, flow accumulations and aspects
  #--------------------------------------------------------------------------------------
  dist_hill <- data.frame('dist2river'=raster::extract(dist2river, hill[, c(1,2)]), 'x'=hill[,1], 'y'=hill[,2])
  elev_hill <- data.frame('elev2river'=raster::extract(elev2river, hill[, c(1,2)]), 'x'=hill[,1], 'y'=hill[,2])
  flow_hill <- data.frame('accum'=raster::extract(accum, hill[, c(1,2)]), 'x'=hill[1], 'y'=hill[2])
  asp_hill <- data.frame('aspect'=raster::extract(aspect, hill[, c(1,2)]), 'x'=hill[1], 'y'=hill[2])
  
  #--------------------------------------------------------------------------------------
  # start calculation
  #--------------------------------------------------------------------------------------
  # reduce precision of dist2river map 
  round_value <- res(dist2river)[1]
  dist_hill <- round(dist_hill/round_value)*round_value
  
  # count all values and there appearance
  ob_mean_catena <- table(dist_hill$dist2river)
  
  # generate representive hillslope data frame Edit - NA values are Ignored (Ashish)
  rep_hill <- data.frame('mean_dist' = as.numeric(names(ob_mean_catena)))
  rep_hill$x <- sapply(names(ob_mean_catena), function(i){mean(elev_hill[dist_hill$dist2river == i,]$x,na.rm=TRUE)})
  rep_hill$y <- as.numeric(sapply(names(ob_mean_catena), function(i){mean(elev_hill[dist_hill$dist2river == i,]$y,na.rm=TRUE)}))
  rep_hill$mean_elev <- as.numeric(sapply(names(ob_mean_catena), function(i){sum(elev_hill$elev2river[dist_hill$dist2river == i] *
                                                                                 sqrt(flow_hill$accum[dist_hill$dist2river == i]),na.rm=TRUE)/
                                                                                 sum(sqrt(flow_hill$accum[dist_hill$dist2river == i]),na.rm=TRUE)}))
  rep_hill$width <- as.numeric(ob_mean_catena)*res(dist2river)[1]
  
  
  ###
 # if geology data is available extract informations - Placeholder for future use
  geo <- li_spatial$geology
  if(exists('geo') & !is.null(geo))
  {
    geo_hill <- data.frame(geo=extract(geo, hill[, c(1,2)]), 'x'=hill[,1], 'y'=hill[,2])
    rep_hill$geo <- as.numeric(sapply(names(ob_mean_catena),
                                      function(i){weighted_flow_accum <- as.integer((flow_hill$accum[dist_hill$dist2river == i] / sum(flow_hill$accum[dist_hill$dist2river == i]))* sum(flow_hill$accum[dist_hill$dist2river == i]))
                                                  names(which.max(table(rep(geo_hill$geo[dist_hill$dist2river == i], weighted_flow_accum))))
                                      }
    ))
    
  }else{rep_hill$geo<-rep(NA, length(rep_hill$mean_elev))}

# if soil data is available extract informations
soil <- li_spatial$soil
if(exists('soil') & !is.null(soil))
{
  soil_hill <- data.frame(soil=raster::extract(soil, hill[, c(1,2)]), 'x'=hill[,1], 'y'=hill[,2])
  rep_hill$soil <- as.numeric(sapply(names(ob_mean_catena),
                                      function(i){
                                        soil_at_dist_indices <- which(dist_hill$dist2river == as.numeric(i)) # Use the original dist2river here

                                        if(length(soil_at_dist_indices) > 0) {
                                          soil_values <- soil_hill$soil[soil_at_dist_indices]
                                          accum_values <- flow_hill$accum[soil_at_dist_indices]

                                          valid_indices <- which(!is.na(soil_values) & !is.na(accum_values) & is.numeric(accum_values) & accum_values >= 0)

                                          if(length(valid_indices) > 0) {
                                            weighted_flow_accum_sum <- sum(accum_values[valid_indices], na.rm = TRUE)
                                            if(weighted_flow_accum_sum > 0) {
                                              weighted_flow_accum <- as.integer((accum_values[valid_indices] / weighted_flow_accum_sum) * weighted_flow_accum_sum)
                                              dom_soil_table <- table(rep(soil_values[valid_indices], weighted_flow_accum))
                                              if(length(dom_soil_table) > 0) {
                                                dominant_soil <- as.numeric(names(which.max(dom_soil_table)))
                                                return(dominant_soil)
                                              } else {
                                                return(NA) # No valid soil after weighting
                                              }
                                            } else {
                                              return(NA) # Sum of weighted accum is zero
                                            }
                                          } else {
                                            return(NA) # No valid soil or accumulation at this distance
                                          }
                                        } else {
                                          return(NA) # No cells at this distance
                                        }
                                      }
  ))

}else{rep_hill$soil <-rep(NA, length(rep_hill$mean_elev))}

# if landuse data is available extract informations
landuse <- li_spatial$landuse
if(exists('landuse') & !is.null(landuse))
{
  landuse_hill <- data.frame(landuse=raster::extract(landuse, hill[, c(1,2)]), 'x'=hill[,1], 'y'=hill[,2])
  rep_hill$landuse <- as.numeric(sapply(names(ob_mean_catena),
                                      function(i){
                                        landuse_at_dist_indices <- which(dist_hill$dist2river == as.numeric(i))

                                        if(length(landuse_at_dist_indices) > 0) {
                                          landuse_values <- landuse_hill$landuse[landuse_at_dist_indices]
                                          accum_values <- flow_hill$accum[landuse_at_dist_indices]

                                          valid_indices <- which(!is.na(landuse_values) & !is.na(accum_values) & is.numeric(accum_values) & accum_values >= 0)

                                          if(length(valid_indices) > 0) {
                                            weighted_flow_accum_sum <- sum(accum_values[valid_indices], na.rm = TRUE)
                                            if(weighted_flow_accum_sum > 0) {
                                              weighted_flow_accum <- as.integer((accum_values[valid_indices] / weighted_flow_accum_sum) * weighted_flow_accum_sum)
                                              dom_landuse_table <- table(rep(landuse_values[valid_indices], weighted_flow_accum))
                                              if(length(dom_landuse_table) > 0) {
                                                dominant_landuse <- as.numeric(names(which.max(dom_landuse_table)))
                                                return(dominant_landuse)
                                              } else {
                                                return(NA) # No valid landuse after weighting
                                              }
                                            } else {
                                              return(NA) # Sum of weighted accum is zero
                                            }
                                          } else {
                                            return(NA) # No valid landuse or accumulation at this distance
                                          }
                                        } else {
                                          return(NA) # No cells at this distance
                                        }
                                      }
  ))

}else{rep_hill$landuse <-rep(NA, length(rep_hill$mean_elev))}
  ###
  # if landuse is available extract informations

  # Distances you’re evaluating (you already use these names)
dists <- as.numeric(names(ob_mean_catena))

# Helper to compute weighted class shares for one distance
dist_list <- lapply(dists, function(d) {
  idx <- which(dist_hill$dist2river == d)
  if (length(idx) == 0) return(NULL)

  lu   <- landuse_hill$landuse[idx]
  acc  <- flow_hill$accum[idx]

  valid <- which(!is.na(lu) & !is.na(acc) & is.finite(acc) & acc >= 0)
  if (length(valid) == 0) return(NULL)

  wsum <- sum(acc[valid])
  if (wsum <= 0) return(NULL)

  # weighted sum per land-use class
  w_per_class <- tapply(acc[valid], lu[valid], sum, default = 0)

  df <- data.frame(
    distance = d,
    class    = as.integer(names(w_per_class)),
    weight   = as.numeric(w_per_class),
    stringsAsFactors = FALSE
  )

  # percent share
  df$pct <- 100 * df$weight / sum(df$weight)

  # dominant class at this distance
  dom_class <- df$class[which.max(df$weight)]
  df$dominant <- (df$class == dom_class)

  df
})

lu_dist <- bind_rows(dist_list)
if (nrow(lu_dist) == 0) stop("No valid land-use distributions could be computed.")

# Optional: provide readable labels for classes (replace with your LUT if you have it)
# class_lut <- c(`1`="Arable", `2`="Forest", `3`="Grass", `4`="Urban", `5`="Water")
# lu_dist$class_label <- recode(as.character(lu_dist$class), !!!class_lut)
lu_dist$class_label <- factor(lu_dist$class)  # numeric classes as factors by default

# Data for dominant annotations (one row per distance)
dom_anno <- lu_dist %>%
  filter(dominant) %>%
  group_by(distance) %>%
  # Place the star slightly above 100% so it’s always visible
  summarise(class_label = first(class_label),
            y = 103, .groups = "drop")

# ---- Plot: one stacked bar per distance, showing distribution across classes ----
jpeg("/out/plots/landuse_distribution.jpg")
gg <- ggplot(lu_dist, aes(x = factor(distance), y = pct, fill = class_label)) +
  geom_col(width = 0.85, color = "grey25") +
  # Star + label for dominant class
  geom_text(data = dom_anno,
            aes(x = factor(distance), y = y, label = paste0("★ ", class_label)),
            fontface = "bold", size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)), limits = c(0, 106),
                     breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%")) +
  labs(
    x = "Distance to river (units of dist_hill$dist2river)",
    y = "Land-use share (%)",
    fill = "Land-use class",
    title = "Distribution of land-use classes along representative hillslope",
    subtitle = "Stacked bars are flow-accumulation–weighted; star marks dominant class per distance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title   = element_text(face = "bold")
  )
print(gg)
dev.off()

jpeg("/out/plots/landuse_distribution2.jpg")
lu_dist$dist_bin <- cut(lu_dist$distance, breaks = pretty(lu_dist$distance, 15), include.lowest = TRUE)

gg <- ggplot(lu_dist, aes(x = dist_bin, y = pct, fill = class_label)) +
  geom_col(width = 0.9, color = "grey25") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Distance to river (binned)", y = "Land-use share (%)", fill = "Class") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(gg)
dev.off()

jpeg("/out/plots/landuse_distribution3.jpg")
ggplot(lu_dist, aes(x = distance, y = pct, fill = class_label)) +
  geom_area(position = "fill", color = NA, alpha = 0.9) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Distance to river", y = "Relative land-use share",
       fill = "Land-use class",
       title = "Land-use composition along hillslope") +
  theme_minimal(base_size = 13)

dev.off()
  # possible to add more information
  
  ###
  #east
  min_east <- rep_hill$x[rep_hill$mean_dist==min(rep_hill$mean_dist)]
  max_east <- rep_hill$x[rep_hill$mean_dist==max(rep_hill$mean_dist)]
  #north
  min_north <- rep_hill$y[rep_hill$mean_dist==min(rep_hill$mean_dist)]
  max_north <- rep_hill$y[rep_hill$mean_dist==max(rep_hill$mean_dist)]
  
  ###
  # calc. mean hillslope aspect
  mean_aspect <- mean(na.omit(asp_hill$aspect))
  # calc. easting and northing
  east <- abs(cos(mean_aspect)*rep_hill$mean_dist)
  north <- abs(sin(mean_aspect)*rep_hill$mean_dist)
  
  ifelse(min_east >= max_east, x <- (min_east-east), x <- (min_east+east))
  ifelse(min_north >= max_north, y <- (min_north-north), y <- (min_north+north))
  
  ###
  # find connected river segment
  xy <- matrix(c(min_east, min_north),1,2)
  river_cell <- raster::extract(stream_id, xy, buffer=100)
  stream_hill_nr <- as.numeric(names(which.max(table(river_cell))))
  
  if(length(rep_hill$mean_dist) > min_dist | area > min_area)
  {
    ###
    #fit a spline to 2d catena and hillslope width
    sp_elev <- smooth.spline(rep_hill$mean_dist, rep_hill$mean_elev, df=length(rep_hill$mean_dist)/freedom)
    sm_spline <- smooth.spline(rep_hill$mean_dist, rep_hill$width, df=length(rep_hill$mean_dist)/freedom)
    
    ###  Checking for short hillslope based on no_rf
    for(id in 1:length(sm_spline$x))
    {
       if(sm_spline$y[id] < quantile(rep_hill$width, no_rf) & rep_hill$mean_dist[id] > quantile(rep_hill$mean_dist, 0.5))
          {idx <- id 
           break
          } else {idx <- id}
    }
    
    short_dist <- rep_hill$mean_dist[1:idx]
    short_elev <- ifelse(predict(sp_elev, rep_hill$mean_dist)$y[1:idx]>0, predict(sp_elev, rep_hill$mean_dist)$y[1:idx], 0)
    short_width <- abs(predict(sm_spline, rep_hill$mean_dist)$y)[1:idx]
    
    ###
    # area correction for short hillslope width
    rep_area <- sum(c(diff(short_dist, lag=1), 0)* short_width)
    
    # area correction by multi. factor area/extent_trapeze
    short_width_area_corr <- short_width * area/rep_area
        
    short_rep_hill <- data.frame('short_dist' = short_dist,
                                 'short_elev' = short_elev,
                                 'short_width' = short_width,
                                 'short_width_corr' = short_width_area_corr,
                                 'east' = x[1:idx],
                                 'north' = y[1:idx],
                                 'geo' = rep_hill$geo[1:idx],
                                 'landuse' = rep_hill$landuse[1:idx],
                                 'soil' = rep_hill$soil[1:idx])
    
    ###
    # plot 2d catena
    if(plot_2d_catena==TRUE)
      {
      plot_2d_profile(rep_hill, sp_elev, short_dist, short_elev, hillslope_nr)
      }
    
    ###
    # plot hillslope width
    if(plot_hillslope_width==TRUE)
      {
      plot_hillslope_width(rep_hill, short_dist, short_width, short_width_area_corr, sm_spline)
      }
     
    li_hill <- list('final_hill' = data.frame('east'=x, 'north'=y,
                                              'dist2river'= rep_hill$mean_dist,
                                              'org_elev' = rep_hill$mean_elev,
                                              'elev2river'= ifelse(predict(sp_elev, rep_hill$mean_dist)$y > 0,
                                                                   predict(sp_elev, rep_hill$mean_dist)$y, 0),
                                              'hillwidth' = rep_hill$width),
                    'all_dist'=dist_hill,
                    'all_elev'=elev_hill,
                    'area'=area, 'stream_id' = stream_hill_nr,
                    'mean_aspect'=mean_aspect,
                    'short_rep_hill'=short_rep_hill,
                    'hillslope_nr'=hillslope_nr
                    )
    cat('  ____________ 
< Hillslope finished successfully!! ')
    
    
  } else {li_hill <- list('final_hill' = data.frame('east'=x, 'north'=y,
                                               'dist2river'= rep_hill$mean_dist,
                                               'org_elev' = rep_hill$mean_elev,
                                               'elev2river'= NA,
                                               'hillwidth' = NA),
                          'all_dist'=dist_hill, 'all_elev'=elev_hill,
                          'area'=area, 'stream_id' = stream_hill_nr,
                          'mean_aspect'=mean_aspect,
                          'hillslope_nr'=hillslope_nr
                         )
  
         print('Hillslope too small or too short!')
         }
  return(li_hill)
}


