#--------------------------------------------------------------------------------------
# This file contains a function based on the hillslope wizards developed
# by Ralf Loritz (2014). The function is used to infer a representative hillslope from
# provided .tif files.

# load Catflow package
library(Catflow)

make_geometry_representative_hillslope <- function(params,data_paths) {
    # loader raster package to open .tif files
    library("raster")

    # load hillslope tool
    source("hillslope_method.R")

    #import spatial data
    flow_accum <- raster(data_paths$flow_accumulation)  # flowaccumulation Attention should be in log scale (Edit: Not working with log - Ashish)
    hillslopes <- raster(data_paths$hillslopes)              # !!! entire basin as one hillslope 
    elev_2_river <- raster(data_paths$elev2river)  # elevation to !!!elev2riv_mod !!! (because of failed calculation areas, gasp are filled with values of SAGA calculation)
    dist_2_river <- raster(data_paths$dist2river)   # distance to river
    dem <- raster(data_paths$filled_dem)                       # digital elevation modell
    aspect <- raster(data_paths$aspect)                 # aspect
    river_id <- raster(data_paths$river_id)             # stream link id

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

    dev.off()

    # project path CATFLOW/in/
    project.path <- "/out/CATFLOW/in/hillgeo"
    system("mkdir -p /out/CATFLOW/in/")
    system("mkdir -p  /out/CATFLOW/in/hillgeo")
    system("chmod 777 /out/CATFLOW")
    system("chmod 777 /out/CATFLOW/in")
    system("chmod 777  /out/CATFLOW/in/hillgeo")

    # transform hillslope to point table
    hillslope_data_frame <- rasterToPoints(hillslopes)

    # create a list of the input maps
    li_spatial <- list("accum" = flow_accum, "hillslopes" = hillslopes, "dem" = dem, "dist2river" = dist_2_river,
                       "elev2river" = elev_2_river, "aspect" = aspect, "hillslope_table" = hillslope_data_frame,
                       "stream_id" = river_id)

    # Select a hillslope from the hillslope raster map (integer value) - halfbasins file
    hillslope_nr <- params$hillslope_id

    # Run hillslope function. no_rf= Percentage of no flow area for region with almost no slope
    # min_dist = minimum number of cells within a hillslope; freedom= freedom of spline function
    pdf("/out/plots/hillslope_plots.pdf")

    hill <- hillslope_tool(hillslope_nr, li_spatial, plot_hillslope_width = TRUE, plot_2d_catena = TRUE,
                           no_rf = params$no_flow_area, min_dist = params$min_cells, min_area = 10000, freedom = 10)

    dev.off()

    #create slope.list for Catflow function make.geometry
    # turn around hillslope catena left to right
    hill$short_rep_hill$east <- rev(hill$short_rep_hill$east)                            # East
    hill$short_rep_hill$north <- rev(hill$short_rep_hill$north)                          # North
    hill$short_rep_hill$short_elev <- rev(hill$short_rep_hill$short_elev + 6)            # Elev + 4 or 6?

    # rectangle, get hillslope width w : A/length = w
    w <- hill$area / hill$short_rep_hill$short_dist[length(hill$short_rep_hill$short_dist)]

    # assign homogenous width to hill object
    hill$short_rep_hill$short_width_corr <- rep(c(w), length(hill$short_rep_hill$short_dist))    # Width

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

    # make geometry from hillslope parameters
    pdf("/out/plots/geometry.pdf")
    out.geom <- make.geometry(topo, make.output = TRUE, project.path = project.path)
    dev.off()

    # save output of make.geometry() for use in other tools
    saveRDS(out.geom, file = "/out/geom.Rds")

    # return to use in workflows
    return(out.geom)
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

write_multipliers <- function(params,data_paths) {
# Output folder
system("mkdir -p /out/CATFLOW/in/")
system("mkdir -p  /out/CATFLOW/in/soil")
system("chmod 777 /out/CATFLOW")
system("chmod 777 /out/CATFLOW/in")
system("chmod 777  /out/CATFLOW/in/soil")

# Assume that rep_hill.geo is already created 
geometry = read.geofile(data_paths$geometry) # This may have to be changed in final version!!!

# write ksmult and thsmult
    write.facmat(
        output.file = "/out/CATFLOW/in/soil/ksmult0.dat", et=geometry$eta, xs=geometry$xsi,
        fac = params$fac_kst
    )
    write.facmat(
        output.file = "/out/CATFLOW/in/soil/thsmult0.dat", et=geometry$eta, xs=geometry$xsi,
        fac = params$fac_ths
    )

}