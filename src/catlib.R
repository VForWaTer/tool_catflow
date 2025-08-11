#--------------------------------------------------------------------------------------
# This file contains a function based on the hillslope wizards developed
# by Ralf Loritz (2014). The function is used to infer a representative hillslope from
# provided .tif files.

# Load required packages
library(Catflow)
library('ggplot2')
library('WVPlots')
library("raster")
library('sf')

source('CustomPlottingFunctions.R')
source("hillslope_method.R")

make_geometry_representative_hillslope <- function(params,data_paths) {
  
    #import spatial data
    flow_accum <- raster(data_paths$flow_accumulation)  # flowaccumulation Attention should be in log scale (Edit: Not working with log - Ashish)
    hillslopes <- raster(data_paths$hillslopes)              # !!! entire basin as one hillslope 
    elev_2_river <- raster(data_paths$elev2river)  # elevation to !!!elev2riv_mod !!! (because of failed calculation areas, gasp are filled with values of SAGA calculation)
    dist_2_river <- raster(data_paths$dist2river)   # distance to river
    dem <- raster(data_paths$filled_dem)                       # digital elevation modell
    aspect <- raster(data_paths$aspect)                 # aspect
    river_id <- raster(data_paths$river_id)             # stream link id
    soil <- raster(data_paths$soil)                     # New soil raster file
    soil_proj <- projectRaster(soil, crs = crs(elev_2_river), method = "ngb") # Use "ngb" for categorical data

    # Read Geopackage if available
 #   hill_gpkg_df <- NULL
  #  if (!is.null(data_paths$hillslopes_vect) && file.exists(data_paths$hillslopes_vect)) {
  #      hillslope_gpkg <- st_read(data_paths$hillslopes_vect, quiet = TRUE)
  #      hill_gpkg_df <- as.data.frame(hillslope_gpkg)
   # }
    
    # plot spatial data
    system("mkdir /out/plots")
    system("chmod 777 /out/plots")

    pdf("/out/plots/spatial_data.pdf")

    plot(hillslopes, main = "Basin")
    plot(river_id, main = "River Network")
    plot(aspect, main = "Aspect")
    plot(flow_accum, main = "Flow Accumulation")
    plot(dem, main = "DEM")
    plot(elev_2_river, main = "Elevation to River")
    plot(dist_2_river, main = "Distance to River")
    plot(soil, main = "Soil Properties")  # Plot the new soil data
    plot(soil_proj)
 
    dev.off()

    # project path CATFLOW/in/
    project.path <- "/out/CATFLOW/in/hillgeo"
    system("mkdir -p /out/CATFLOW/in/")
    system("mkdir -p  /out/CATFLOW/in/hillgeo")
    system("mkdir -p  /out/CATFLOW/in/soil") 
    system("chmod 777 /out/CATFLOW")
    system("chmod 777 /out/CATFLOW/in")
    system("chmod 777  /out/CATFLOW/in/hillgeo")
    system("chmod 777  /out/CATFLOW/in/soil")

    # transform hillslope to point table
    hillslope_data_frame <- rasterToPoints(hillslopes)

    # Remove decimal places (performance)
    elev_2_river <- floor(elev_2_river)
    dist_2_river <- floor(dist_2_river)
    flow_accum <-  floor(flow_accum)

    # create a list of the input maps
    li_spatial <- list("accum" = flow_accum, "hillslopes" = hillslopes, "dem" = dem, "dist2river" = dist_2_river,
                       "elev2river" = elev_2_river, "aspect" = aspect, "hillslope_table" = hillslope_data_frame,
                       "stream_id" = river_id, "soil" = soil_proj)  # Add soil to the spatial data list

    # Select a hillslope from the hillslope raster map (integer value) - halfbasins file
    hillslope_nr <- params$hillslope_id

    # Run hillslope function. no_rf= Percentage of no flow area for region with almost no slope
    # min_dist = minimum number of cells within a hillslope; freedom= freedom of spline function
    pdf("/out/plots/hillslope_plots.pdf")

    hill <- hillslope_tool(
        hillslope_nr, li_spatial, plot_hillslope_width = TRUE, plot_2d_catena = TRUE,
        no_rf = params$no_flow_area, min_dist = params$min_cells,
        min_area = params$min_area, freedom = params$freedom
    )

    dev.off()

    if (!is.null(hill)) {
        pdf("/out/plots/energy_distribution.pdf")
        Data <- data.frame(Distance = hill$all_dist$dist2river, Elevation = hill$all_elev$elev2river)
        ScatterHist(
          Data, "Distance", "Elevation", title = "Energy distribution plot", hill = hill,
          smoothmethod = "none", contour = FALSE, point_color = "#006d2c",
          hist_color = "#6baed6", density_color = "red", point_alpha = 0.05,
          minimal_labels = TRUE
        ) + theme_bw()
        dev.off()
        # Also save as JPEG
        jpeg("/out/plots/energy_distribution.jpg")
        ScatterHist(
          Data, "Distance", "Elevation", title = "Energy distribution plot", hill = hill,
          smoothmethod = "none", contour = FALSE, point_color = "#006d2c",
          hist_color = "#6baed6", density_color = "red", point_alpha = 0.05,
          minimal_labels = TRUE
        ) + theme_bw()
        dev.off()
        # Also save as JPEG with higher quality and resolution
        jpeg("/out/plots/energy_distribution_hq.jpg", width = 4000, height = 3000, res = 600, quality = 100)
        ScatterHist(
          Data, "Distance", "Elevation", title = "Energy distribution plot", hill = hill,
          smoothmethod = "none", contour = FALSE, point_color = "#006d2c",
          hist_color = "#6baed6", density_color = "red", point_alpha = 0.05,
          minimal_labels = TRUE
        ) + theme_bw()
        dev.off()

        #create slope.list for Catflow function make.geometry
        # turn around hillslope catena left to right
        hill$short_rep_hill$east <- rev(hill$short_rep_hill$east)                            # East
        hill$short_rep_hill$north <- rev(hill$short_rep_hill$north)                          # North
        hill$short_rep_hill$short_elev <- rev(hill$short_rep_hill$short_elev + 6)            # Elev + 4 or 6?

        # rectangle, get hillslope width w : A/length = w
        w <- hill$area / hill$short_rep_hill$short_dist[length(hill$short_rep_hill$short_dist)]

        # assign homogenous width to hill object if constant_width parameter is true
        if (params$constant_width) {
            hill$short_rep_hill$short_width_corr <- rep(c(w), length(hill$short_rep_hill$short_dist))    # Width
        }

        # Assign values to x based on the value of params$hill_type
        if (params$hill_type == "constant") {
        htyp <- 1
        } else if (params$hill_type == "cake") {
        htyp <- 2
        } else if (params$hill_type == "variable") {
        htyp <- 3
        } 
        
        # slope.list as input for make_geometry
        topo <- list(
            xh = hill$short_rep_hill$east,              # East
            yh = hill$short_rep_hill$north,             # North
            zh = hill$short_rep_hill$short_elev,        # Elev
            bh = hill$short_rep_hill$short_width_corr,  # Width
            dist = hill$short_rep_hill$short_dist,
            tot.area = hill$area,                       # Area
            htyp = htyp,                                   # Hillslope type
            dyy = params$depth,                                    # Thickness of profile [m] average of drillings 2.1 m
            xsi = seq(0, 1, length = max(hill$short_rep_hill$short_dist) + 1),  # in order to get every 1m a node length of hillslope + 1             # discretitation
            eta =  c(seq(0,0.625,length=6),seq(0.7,0.85, length=3),seq(0.875,1, length=6)),       # eta starts at the bottom, upper 50 cm, dx=10cm, 50-400cm dx=25 cm
            out.file = "rep_hill.geo"           # outpath
        )

      #  if (!is.null(hill_gpkg_df) && params$hillslope_id != -1) {
       #     precid_val <- hill_gpkg_df$hillslope_new_id[hill_gpkg_df$hillslope_id == params$hillslope_id][1]
      #      if (is.na(precid_val)) precid_val <- 1
       # } else {
      #      precid_val <- 1
       # }
        # make geometry from hillslope parameters
        precid_val <- 1
        out.geom <- make.geometry(topo,numh =precid_val, make.output = TRUE, project.path = project.path)

        pdf("/out/plots/geometry.pdf")
        plot.catf.grid(out.geom$sko, out.geom$hko, val=out.geom$hko,boundcol = 1)
        dev.off()

        # 0. Check if hill$short_rep_hill has at least one non-NA value for soil
        if (any(!is.na(hill$short_rep_hill$soil))) {

            # Extract the unique soil types from hill$short_rep_hill
            unique_soil_types <- unique(hill$short_rep_hill$soil)

            # Initialize the file content with the first line
            file_content <- paste(length(unique_soil_types), "0", sep = " ")

            # Initialize a variable to track the starting xsi for the next soil type
            current_xsi_start <- 0.0

            # Loop through each unique soil type
            for (soil_type in unique_soil_types) {
                # Find the minimum and maximum 'short_dist' values for the current soil type
                short_dist_values <- hill$short_rep_hill$short_dist[hill$short_rep_hill$soil == soil_type]
                short_dist_min <- min(short_dist_values)
                short_dist_max <- max(short_dist_values)

                # Calculate the overall horizontal extent
                short_dist_min_overall <- min(hill$short_rep_hill$short_dist)
                short_dist_max_overall <- max(hill$short_rep_hill$short_dist)

                # Calculate relative xsi coordinates
                relative_xsi_min <- current_xsi_start
                relative_xsi_max <- (short_dist_max - short_dist_min_overall) / (short_dist_max_overall - short_dist_min_overall)

                # Append the line to the file content string
                file_content <- paste0(
                    file_content, "\n",
                    "0.0 1.0 ", sprintf("%.1f", relative_xsi_min), " ", sprintf("%.1f", relative_xsi_max), " ", soil_type
                )

                # Update the starting xsi for the next soil type
                current_xsi_start <- relative_xsi_max
            }

            # Write the content to a file
            output_dir <- "/out/CATFLOW/in/soil"
            output_file <- file.path(output_dir, "soil.dat")
            write(file_content, file = output_file)

        } else {
            cat("hill$short_rep_hill$soil contains only NA values. No soil definition file will be created.\n")
        }
  
        # Parse file_content and convert it into a matrix
        if (exists("file_content")) {
            # Split the file_content into lines
            lines <- strsplit(file_content, "\n")[[1]]
            
            # Skip the first line (header) and process the remaining lines
            soil_data <- do.call(rbind, lapply(lines[-1], function(line) {
                values <- strsplit(line, " ")[[1]]
                as.numeric(values)
            }))
            
            # Align the soil matrix with the dimensions of out.geom$sko
            sko_dims <- dim(out.geom$sko)  # Get dimensions of sko
            soil_matrix <- matrix(NA, nrow = sko_dims[1], ncol = sko_dims[2])  # Initialize matrix
            
            # Fill the matrix with soil type values
            for (i in 1:nrow(soil_data)) {
                # Map vertical and horizontal grid indices
                row_start <- round(soil_data[i, 1] * (sko_dims[1] - 1)) + 1
                row_end <- round(soil_data[i, 2] * (sko_dims[1] - 1)) + 1
                col_start <- round(soil_data[i, 3] * (sko_dims[2] - 1)) + 1
                col_end <- round(soil_data[i, 4] * (sko_dims[2] - 1)) + 1

                # Ensure indices are within bounds of the domain
                row_start <- max(1, min(row_start, sko_dims[1]))
                row_end <- max(1, min(row_end, sko_dims[1]))
                col_start <- max(1, min(col_start, sko_dims[2]))
                col_end <- max(1, min(col_end, sko_dims[2]))

                # Assign the soil type to all rows and columns in the range
                soil_matrix[row_start:row_end, col_start:col_end] <- soil_data[i, 5]
            }
            
            # Fill NA values in the matrix with a default soil type (e.g., 0 for no data)
            soil_matrix[is.na(soil_matrix)] <- 0
        }

        # save output of make.geometry() for use in other tools
        saveRDS(out.geom, file = "/out/geom.Rds")
        
        pdf("/out/plots/soil_types.pdf")
        plot.catf.grid(out.geom$sko, out.geom$hko, val=soil_matrix)
        dev.off()

        # write mulipliers file
        system("mkdir -p /out/CATFLOW/in/")
        system("mkdir -p  /out/CATFLOW/in/soil")
        system("chmod 777 /out/CATFLOW")
        system("chmod 777 /out/CATFLOW/in")
        system("chmod 777  /out/CATFLOW/in/soil")

        #defaults to 1 change manually if detailed soil data is available   
        fact_mult <- 1
        write.facmat(
        output.file = "/out/CATFLOW/in/soil/ksmult0.dat", et=out.geom$eta, xs=out.geom$xsi,
        fac = fact_mult
        )
        write.facmat(
            output.file = "/out/CATFLOW/in/soil/thsmult0.dat", et=out.geom$eta, xs=out.geom$xsi,
            fac = fact_mult
        )

        # write surface nodes file
        system("mkdir -p  /out/CATFLOW/in/landuse")
        system("chmod 777  /out/CATFLOW/in/landuse")
        
        write.surface.pob(output.file = "/out/CATFLOW/in/landuse/surface.pob", 
                  xs = out.geom$xsi, lu = 1, precid = precid_val, climid = 1, 
                  windid = rep(1, 4))

        # return to use in workflows
        return(out.geom)

    } else {
        stop("Hillslope tool did not return a valid hill object.")
    }
}
