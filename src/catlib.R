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
    landuse <- raster(data_paths$landuse)                 # landuse

    # Read Geopackage
    hillslope_gpkg <- st_read(data_paths$hillslopes_vect)
    hill_gpkg_df = as.data.frame(hillslope_gpkg)
    
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
    plot(landuse, main = "Land Use")  # Plot the land use data
 
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
                       "stream_id" = river_id, "soil" = soil_proj, "landuse" = landuse)  # Add soil and landuse to the spatial data list

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

    short_rep_hill <- hill$short_rep_hill
    # save as RDS (and optionally as RData)
    save(short_rep_hill, file = file.path("/out/CATFLOW/in/hillgeo", "hill.RData"))

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
        hill$short_rep_hill$short_elev <- rev(hill$short_rep_hill$short_elev)            #Elev
        hill$short_rep_hill$short_width <- rev(hill$short_rep_hill$short_width)            #Width
        hill$short_rep_hill$short_width_corr <- rev(hill$short_rep_hill$short_width_corr)            #Width

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

        # This is a workaround to ensure that hillslope id is assigned sequentially as CATFLOW doesnt explictily allow the definition of a hillslope id 
        if (params$hillslope_id != -1) {
            precid_val <- hill_gpkg_df$hillslope_new_id[hill_gpkg_df$hillslope_id == params$hillslope_id][1]
        } else {
            precid_val <- 1
        }
        # make geometry from hillslope parameters
        out.geom <- make.geometry(topo,numh =precid_val, make.output = TRUE, project.path = project.path)

        pdf("/out/plots/geometry.pdf")
        plot.catf.grid(out.geom$sko, out.geom$hko, val=out.geom$hko,boundcol = 1)
        dev.off()

# --- Inputs assumed: out.geom$xsi in [0,1], hill$short_rep_hill with short_dist and soil ---

        # 0) Quick checks
        stopifnot(is.numeric(out.geom$xsi), length(out.geom$xsi) > 0)
        hrh <- hill$short_rep_hill
        if (!("short_dist" %in% names(hrh)) || !("soil" %in% names(hrh))) {
        stop("hill$short_rep_hill must contain 'short_dist' and 'soil' columns.")
        }
        if (!any(!is.na(hrh$soil))) {
        stop("No non-NA soils found in hill$short_rep_hill$soil.")
        }

        # 1) Relative position of segments (0..1)
        x_rel_seg <- hrh$short_dist / max(hrh$short_dist, na.rm = TRUE)

        # 2) Ensure surface nodes (xsi) are sorted and remember original order (not strictly needed for files)
        ord <- order(out.geom$xsi)
        xsi_sorted <- out.geom$xsi[ord]

        # 3) Nearest segment per node (fast, no loops)
        iL <- pmax(1L, pmin(findInterval(xsi_sorted, x_rel_seg), length(x_rel_seg) - 1L))
        iR <- iL + 1L
        choose_right <- abs(xsi_sorted - x_rel_seg[iR]) < abs(xsi_sorted - x_rel_seg[iL])
        nearest_idx <- ifelse(choose_right, iR, iL)

        soil_vec_sorted <- hrh$soil[nearest_idx]

        # Handle possible NA soils by simple forward/back fill
        if (anyNA(soil_vec_sorted)) {
        # forward fill then backward fill
        for (k in seq_along(soil_vec_sorted)) {
            if (is.na(soil_vec_sorted[k]) && k > 1) soil_vec_sorted[k] <- soil_vec_sorted[k-1]
        }
        for (k in length(soil_vec_sorted):1) {
            if (is.na(soil_vec_sorted[k]) && k < length(soil_vec_sorted)) soil_vec_sorted[k] <- soil_vec_sorted[k+1]
        }
        # if still NA (all were NA), stop
        if (anyNA(soil_vec_sorted)) stop("Soil mapping produced only NA values after fill.")
        }

        # 4) Build node "cells" (intervals) via midpoints between xsi nodes
        #    bounds has length N+1: [0, midpoint(1-2), ..., midpoint(N-1 - N), 1]
        N <- length(xsi_sorted)
        if (N == 1L) {
        bounds <- c(0.0, 1.0)
        } else {
        mids <- head(xsi_sorted, -1) + diff(xsi_sorted) / 2
        bounds <- c(0.0, mids, 1.0)
        }

        # 5) Collapse consecutive nodes with the same soil into bands
        r <- rle(soil_vec_sorted)
        run_lengths <- r$lengths
        run_values  <- r$values

        # Start indices of runs in 1..N
        run_starts <- cumsum(c(1, head(run_lengths, -1)))
        run_ends   <- cumsum(run_lengths)

        # Each run spans [bounds[start], bounds[end+1]]
        xsi_min_vec <- bounds[run_starts]
        xsi_max_vec <- bounds[run_ends + 1]
        soil_id_vec <- run_values

        # 6) Compose file content
        #    IMPORTANT: header should list the number of bands (intervals), not unique soils,
        #    because a soil may appear in multiple disjoint bands.
        num_bands <- length(soil_id_vec)
        file_content <- paste(num_bands, "0")  # your existing "N 0" header

        # Append one line per band: "0.0 1.0 <xsi_min> <xsi_max> <soil_id>"
        # Match your formatting (one decimal); adjust to "%.3f" if you want finer granularity.
        lines <- sprintf("0.0 1.0 %.1f %.1f %s",
                        xsi_min_vec, xsi_max_vec, as.character(soil_id_vec))

        file_content <- paste(file_content, paste(lines, collapse = "\n"), sep = "\n")

        # 7) (Optional) write to file
        writeLines(file_content, "/out/CATFLOW/in/soil/soil.dat")  # or your desired path

        # save output of make.geometry() for use in other tools
        saveRDS(out.geom, file = "/out/geom.Rds")

        # write mulipliers file
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
        
        # Find hillslope_new_id from hill_gpkg_df where hillslope_id matches params$hillslope_id
        precid_val <- hill_gpkg_df$hillslope_new_id[hill_gpkg_df$hillslope_id == params$hillslope_id][1]
        hrh <- hill$short_rep_hill
        x_rel_seg <- hrh$short_dist / max(hrh$short_dist)
        # 2) For every surface node, find nearest segment index
        nearest_seg_idx <- sapply(out.geom$xsi, function(x) which.min(abs(x_rel_seg - x)))

        # 3) Take the landuse ID from the segment table
        landuse_vec <- hrh$landuse[nearest_seg_idx]

        landuse_vec = rev(landuse_vec)  # Reverse to match CATFLOW orientation

        write.surface.pob(output.file = "/out/CATFLOW/in/landuse/surface.pob", 
                  xs = out.geom$xsi, lu = landuse_vec, precid = precid_val, climid = 1, 
                  windid = rep(1, 4))

        # return to use in workflows
        return(out.geom)

    } else {
        stop("Hillslope tool did not return a valid hill object.")
    }
}

define_run_printouts <- function(params,data_paths) {
# Output folder
system("mkdir -p /out/CATFLOW/in/")
system("mkdir -p  /out/CATFLOW/in/control")
system("chmod 777 /out/CATFLOW")
system("chmod 777 /out/CATFLOW/in")
system("chmod 777  /out/CATFLOW/in/control")

# Assign values to x based on the value of params$hill_type
    if (params$time.unit == "hourly") {
    typ <- 'h'
    } else if (params$hill_type == "seconds") {
    typ <- 's'
    } 
# write printout times
    write.printout(
        output.file = "/out/CATFLOW/in/control/printout.prt",
        start.time = params$start.time,
        end.time = params$end.time,
        intervall = params$interval,
        time.unit = typ,
        flag = params$flag
    )

}
