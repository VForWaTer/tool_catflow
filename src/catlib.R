#--------------------------------------------------------------------------------------
# This file essentially maps the Catflow preprocessing functions from the
# Catflow-R-Pacakge and makes them available in the container.
# "Real" workflows which combine the different functions can be found in lib.R
# This file additionally contains a function based on the hillslope wizards developed
# by Ralf Loritz (2014). The function is used to infer a representative hillslope from
# provided .tif files.

# load Catflow package
library(Catflow)

catlib_make_geometry <- function(params) {
    # output file is always saved to /out/CATFLOW/in/
    params$out.file <- "geometry.geo"

    # create project folder
    project.path <- "/out/CATFLOW/in"
    system("mkdir -p /out/CATFLOW/in/")
    system("chmod 777 /out/CATFLOW")
    system("chmod 777 /out/CATFLOW/in")

    # create folder for plots
    system("mkdir -p /out/plots")
    system("chmod 777 /out/plots")

    # run function make.geometry() with params as input parameters
    pdf("/out/plots/geometry.pdf")
    output <- make.geometry(params, plottitle = "Model Geometry", project.path = project.path)
    dev.off()

    # save output of make.geometry() for use in other tools
    saveRDS(output, file = "/out/geom.Rds")
}

catlib_write_facmat <- function(params) {
    # load geometry generated by the tool make_geometry and attach
    geom <- readRDS(params$geometry)
    attach(geom)

    # create project folder
    system("mkdir -p /out/CATFLOW/in")
    system("chmod 777 /out/CATFLOW")
    system("chmod 777 /out/CATFLOW/in")

    # multipliers for scaling saturated hydraulic conductivity
    write.facmat(output.file = "/out/CATFLOW/in/ksmult.dat")

    # multipliers for scaling saturated water content / porosity
    write.facmat(output.file = "/out/CATFLOW/in/thsmult.dat")

    # write initial conditions for hydraulic conductivity
    if (params$write_soilhyd_ini) {
       write.facmat(output.file = "/out/CATFLOW/in/soilhyd.ini",
                    headr = paste("PSI ", 0, 1, length(eta), length(xsi), 1))
    }

    # write soiltype identifiers
    if (params$write_soil_types) {
        write.facmat(output.file = "/out/CATFLOW/in/soils.bod",
                     headr = paste("BODEN", length(eta), length(xsi), 1),
                     fac = matrix(c(rep(1, ceiling(length(eta) / 2)),
                                    rep(2, floor(length(eta) / 2))),
                                  nrow = length(eta), ncol = length(xsi)))
    }
}

catlib_write_precip <- function(params) {
    # create project folder
    system("mkdir -p /out/CATFLOW/in")
    system("chmod 777 /out/CATFLOW")
    system("chmod 777 /out/CATFLOW/in")

    # write precipitation data with params as input
    write.precip(
        raindat = params$raindat,
        output.file = "/out/CATFLOW/in/rain.dat",
        start.time = params$start.time,
        time.unit = params$time.unit,
        faktor.p = params$faktor.p
    )

    # create folder for plots
    system("mkdir -p /out/plots")
    system("chmod 777 /out/plots")

    # plot rainfall data
    pdf("/out/plots/raindat.pdf")
    plot(params$raindat, t = "s", xlab = "time", ylab = "precipitation")
    dev.off()
}

catlib_write_climate <- function(params) {
    # create project folder
    system("mkdir -p /out/CATFLOW/in")
    system("chmod 777 /out/CATFLOW")
    system("chmod 777 /out/CATFLOW/in")

    # write climate data with params as input
    write.climate(
        climadat = params$climadat,
        output.file = "/out/CATFLOW/in/clima.dat",
        start.time = params$start.time,
        time.unit = params$time.unit,
        rBilart = params$rBilart,
        ref.height = params$ref.height,
        sw0 = params$sw0,
        sw1 = params$sw1,
        sw2 = params$sw2,
        trueb = params$trueb,
        truebf = params$truebf,
        NA.flag = params$NA.flag
    )
}

catlib_write_printout <- function(params) {
    # create project folder
    system("mkdir -p /out/CATFLOW/in")
    system("chmod 777 /out/CATFLOW")
    system("chmod 777 /out/CATFLOW/in")

    # write printout times
    write.printout(
        output.file = "/out/CATFLOW/in/printout.prt",
        start.time = params$start.time,
        end.time = params$end.time,
        intervall = params$intervall,
        time.unit = params$time.unit,
        flag = params$flag
    )
}

catlib_write_surface_pob <- function(params) {
    # load geometry generated by the tool make_geometry and attach
    geom <- readRDS(params$geometry)
    attach(geom)

    # create project folder
    system("mkdir -p /out/CATFLOW/in")
    system("chmod 777 /out/CATFLOW")
    system("chmod 777 /out/CATFLOW/in")

    # drop geometry from params (unused in write.surface.pob)
    params$geometry <- NULL

    # write surface attributes
    write.surface.pob(
        output.file = "/out/CATFLOW/in/surface.pob",
        xs = xsi,
        lu = params$lu,
        precid = params$precid,
        climid = params$climid,
        windid = params$windid
    )
}

catlib_write_control <- function(params) {
    # create the project folder and set permissions
    system("mkdir -p /out/CATFLOW")
    system("chmod 777 /out/CATFLOW")

    # write project controll file
    write.control(
        output.file = "run_cat.in",
        project.path = "/out/CATFLOW",
        start.date = params$start.time,
        end.date = params$end.time,
        slope.in.list = list(
            slope1 = list(
                geo.file = "geometry.geo",
                soil.file = "soils.bod",
                ks.fac = "ksmult.dat",
                ths.fac = "thsmult.dat",
                macro.file = "profil.mak",
                cv.file = "cont_vol.cv",
                ini.file = "soilhyd.ini",
                print.file = "printout.prt",
                surf.file = "surface.pob",
                bc.file = "boundary.rb"
            )
        )
    )

    # write main control file
    write.CATFLOW.IN(
        control.files = "run_cat.in",
        project.path = "/out/CATFLOW"
    )
}

catlib_complete_file_structure <- function(params) {
    # create the project folder CATFLOW and the in/ folder and set permissions
    system("mkdir -p /out/CATFLOW/in")
    system("chmod 777 /out/CATFLOW")
    system("chmod 777 /out/CATFLOW/in")

    # macro.file: "profil.mak"
    cat(paste("1 0 2", "ari", "0.00 1.00 0.00 1.00 1 1.00 1.00 ", sep = "\n"),  file = "/out/CATFLOW/in/profil.mak")

    # cv.file: "cont_vol.cv"
    cat(paste("1", "0.8 0.9 0.98 1.0", sep = "\n"), file = "/out/CATFLOW/in/cont_vol.cv")

    # bc.file: "boundary.rb"
    cat(paste("L", "1 0", "0. 1. 0", " ",
              "R", "1 0", "0. 1. -10", " ",
              "T", "1 0", "0. 1. -99 ", " ",
              "B", "1 0", "0. 1. 0", " ",
              "S", "1 0", "0. 1. 0. 1. -99", " ",
              "M", 0, sep = "\n"),
              file = "/out/CATFLOW/in/boundary.rb")

    # winddir.file: "winddir.def"
    cat(paste("4", "240 0.81", " 50 0.78", " 80 0.97", "220 0.94", sep = "\n"), file = "/out/CATFLOW/in/winddir.def")

    # soildef.file: soils.def
    # writes a soil type definition for two soil types:
    cat(paste("2", "1 Loamy Sand, porosity 0.55, bulk dens 1 g/cm3",
              "1 800 1. 1. 1e-4 0.5 0.34 0.11 20. 0.70 0.050 1. 1. 1.",
              "4.05e-5 0.55 0.06 12.40 2.28 -6.00 8.00 1000.00 0.80",
              "0. 0. 0.", "0. 0. 0.", "0. 0. 0.",
              "2 Sandy Clay Loam (30% S, 40 % U; 30 % T)",
              "1 800 1. 1. 1e-4 0.5 0.34 0.11 20. 0.70 0.050 1. 1. 1.",
              "3.42e-6 0.48 0.08 0.96 1.5 -6.00 8.00 1200.00 0.80",
              "0. 0. 0.", "0. 0. 0.", "0. 0. 0.", sep = "\n"),
              file = "/out/CATFLOW/in/soils.def")

    # timeser.file: timeser.def
    cat(paste("PREC", "1", "in/rain.dat", "",
              "BC", "0", "", "SINKS", "0", "", "SOLUTE", "0", "",
              "LAND-USE", "in/landuse/lu_ts.dat", "",
              "CLIMATE", "1", "in/clima.dat", "", sep = "\n"),
              file = "/out/CATFLOW/in/timeser.def")

    # file related to landuse specification
    # create landuse subdirectory and set permissions
    system("mkdir -p /out/CATFLOW/in/landuse")
    system("chmod 777 /out/CATFLOW/in/landuse")

    # pointer to land-use parameters
    cat(paste("3", "coniferous forest", "in/landuse/conif.par", sep = "              "), file = "/out/CATFLOW/in/landuse/lu_file.def")

    # time-series of land-use parameters
    cat(paste("01.01.2004 00:00:00.00", "in/landuse/lu_set1.dat", "01.01.2005 00:00:00.00", sep = "\n"), file = "/out/CATFLOW/in/landuse/lu_ts.dat")

    # land-use parameters
    cat(paste(
         paste("10", "KST", "MAK", "BFI", "BBG", "TWU", "PFH",
               "PALB", "RSTMIN", "WP_BFW", "F_BFW", sep = "    "),
         "0.    3.    1.    5.    0.95    5.0    5.0    0.15    1.    1.    1.",
         paste(c(" 1", "366"),
               "2.    1.    1.    1.0    1.0    1.0    1.0    546.    0.05    30.",
         sep = "    ", collapse = "\n"), sep = "\n"),
     file = "/out/CATFLOW/in/landuse/conif.par")

    # pointer to surface node attributes
    cat(paste(1, "33    3    %coniferous forest", sep = "\n"), file = "/out/CATFLOW/in/landuse/lu_set1.dat")
}

catlib_make_geometry_representative_hillslope <- function(params) {
    # loader raster package to open .tif files
    library("raster")

    # load hillslope tool
    source("hillslope_method.R")

    #import spatial data
    flow_accum <- raster(params$flow_accumulation)  # flowaccumulation Attention should be in log scale (Edit: Not working with log - Ashish)
    hillslopes <- raster(params$basin)              # !!! entire basin as one hillslope 
    elev_2_river <- raster(params$elevation2river)  # elevation to !!!elev2riv_mod !!! (because of failed calculation areas, gasp are filled with values of SAGA calculation)
    dist_2_river <- raster(params$distance2river)   # distance to river
    dem <- raster(params$dem)                       # digital elevation modell
    aspect <- raster(params$aspect)                 # aspect
    river_id <- raster(params$river_id)             # stream link id

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
    project.path <- "/out/CATFLOW/in"
    system("mkdir -p /out/CATFLOW/in/")
    system("chmod 777 /out/CATFLOW")
    system("chmod 777 /out/CATFLOW/in")

    # transform hillslope to point table
    hillslope_data_frame <- rasterToPoints(hillslopes)

    # create a list of the input maps
    li_spatial <- list("accum" = flow_accum, "hillslopes" = hillslopes, "dem" = dem, "dist2river" = dist_2_river,
                       "elev2river" = elev_2_river, "aspect" = aspect, "hillslope_table" = hillslope_data_frame,
                       "stream_id" = river_id)

    # Select a hillslope from the hillslope raster map (integer value) - halfbasins file
    hillslope_nr <- 1

    # Run hillslope function. no_rf= Percentage of no flow area for region with almost no slope
    # min_dist = minimum number of cells within a hillslope; freedom= freedom of spline function
    pdf("/out/plots/hillslope_plots.pdf")

    hill <- hillslope_tool(hillslope_nr, li_spatial, plot_hillslope_width = TRUE, plot_2d_catena = TRUE,
                           no_rf = 0.37, min_dist = 10, min_area = 10000, freedom = 10)

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

    # slope.list as input for make_geometry
    topo <- list(
        xh = hill$short_rep_hill$east,              # East
        yh = hill$short_rep_hill$north,             # North
        zh = hill$short_rep_hill$short_elev,        # Elev
        bh = hill$short_rep_hill$short_width_corr,  # Width
        dist = hill$short_rep_hill$short_dist,
        tot.area = hill$area,                       # Area
        htyp = 1,                                   # Hillslope type
        dyy = 1,                                    # Thickness of profile [m] average of drillings 2.1 m
        xsi = seq(0, 1, length = max(hill$short_rep_hill$short_dist) + 1),  # in order to get every 1m a node length of hillslope + 1             # discretitation
        eta = c(seq(0, 0.60, length = 4), seq(0.80, 1, length = 11)),       # eta starts at the bottom, upper 50 cm, dx=10cm, 50-400cm dx=25 cm
        out.file = "geometry.geo"           # outpath
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