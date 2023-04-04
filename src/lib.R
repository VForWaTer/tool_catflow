# source catlib.R
source("catlib.R")


lib_preprocess_catflow <- function(params) {
    ###
    # CATFLOW preprocessing workflow, creates all the files
    # necessary to start a CATFLOW run.
    ###

    # (1) make geometry from .tif files
    catlib_make_geometry_representative_hillslope(params)

    # geom.Rds produced above as input for other tools
    params$geometry <- "/out/geom.Rds"

    # (2) create ksmult.dat, thsmult.dat, soilhyd.ini, soils.bod
    params$write_soilhyd_ini <- TRUE
    params$write_soil_types <- TRUE

    catlib_write_facmat(params)

    # (3) write precipitation data file
    catlib_write_precip(params)

    # (4) write climate data file
    catlib_write_climate(params)

    # (5) write printout file
    catlib_write_printout(params)

    # (6) write surface pob file
    catlib_write_surface_pob(params)

    # (7) write controll file
    catlib_write_control(params)

    # (8) complete the CATFLOW file structure
    catlib_complete_file_structure(params)
}


lib_timeseries_to_catflow_precip <- function(data, start.time, time.unit, faktor.p) {
    ###
    # Convert precipitation time series data into the
    # format required by CATFLOW.
    #
    # Always give the precipitation data (data) in mm/h.
    ###

    # Convert the time column to a datetime object
    data[, 1] <- as.POSIXct(data$time, format = "%Y-%m-%d %H:%M:%S")

    if (start.time) {
        # Create a new column with the hours since the given start.time
        data$hours <- as.numeric(difftime(data[, 1], start.time, units = "hours"))
    } else {
        # Create a new column with the hours since the start of the time series
        data$hours <- as.numeric(difftime(data[, 1], min(data$time), units = "hours"))
    }

    catlib_write_precip(
        raindat = data,
        start.time = start.time,
        time.unit = time.unit,
        faktor.p = faktor.p
    )
}
