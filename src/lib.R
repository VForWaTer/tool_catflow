# source catlib.R
source("catlib.R")

# implement "workflow scripts"
lib_preprocess_catflow <- function(params) {
    # (1) make geometry from .tif files
    catlib_make_geometry_representative_hillslope(params)

    # geom.Rds produced above as input for other ttols
    params$geometry <- "/out/geom.Rds"

    # (2) create ksmult.dat, thsmult.dat, soilhyd.ini, soils.bod
    params$write_soilhyd_ini <- TRUE
    params$write_soil_types <- TRUE

    catlib_write_facmat(params)
}
