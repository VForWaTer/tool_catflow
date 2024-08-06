# load json2aRgs for parameter parsing
library(json2aRgs)

# load Catflow package
library(Catflow)

# get the call parameters for the tool
params <- get_parameters()

# check if a toolname was set in env
toolname <- tolower(Sys.getenv("TOOL_RUN"))

# make functions from Catflow-R-Package available
source("catlib.R")

# workflow functions in lib.R
source("lib.R")

# default tool make geometry -> future decission for a default tool
if (toolname == "") {
    toolname <- "make_geometry"
}

# Switch for the different tools available in this package
if (toolname == "make_geometry") {
    catlib_make_geometry(params)
} else if (toolname == "write_facmat") {
    catlib_write_facmat(params)

} else if (toolname == "write_precip") {
    catlib_write_precip(params)

} else if (toolname == "write_climate") {
    catlib_write_climate(params)

} else if (toolname == "write_printout") {
    catlib_write_printout(params)

} else if (toolname == "write_surface_pob") {
    catlib_write_surface_pob(params)

} else if (toolname == "write_control") {
    catlib_write_control(params)

} else if (toolname == "complete_file_structure") {
    catlib_complete_file_structure(params)

} else if (toolname == "make_geometry_representative_hillslope") {
    catlib_make_geometry_representative_hillslope(params)

} else if (toolname == "preprocess_catflow") {
    lib_preprocess_catflow(params)

} else if (toolname == "timeseries_to_catflow_precip") {
    lib_timeseries_to_catflow_precip(
        data = params$data,
        start.time = params$start.time,
        time.unit = params$time.unit,
        faktor.p = params$faktor.p
        )
        
} else {
    # in any other case, the tool was invalid or not configured
    print(paste("[", Sys.time(), "] Either no TOOL_RUN environment variable available, or '", toolname, "' is not valid.\n", sep = ""))
}