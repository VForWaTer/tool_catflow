# source catlib.R
source("catlib.R")

# implement "workflow scripts"
lib_preprocess_catflow <- function(params) {
    # make geometry from .tif files
    geom <- catlib_make_geometry_representative_hillslope(params)
}
