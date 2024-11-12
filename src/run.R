# load json2aRgs for parameter parsing
library(json2aRgs)

# load Catflow package
library(Catflow)

# get the call parameters for the tool
params <- get_parameters()

#get data paths
data_paths = get_data()

# check if a toolname was set in env
toolname <- tolower(Sys.getenv("TOOL_RUN"))

# make functions from Catflow-R-Package available
source("catlib.R")

# Switch for the different tools available in this package
if (toolname == "make_representative_hillslope") {
    make_geometry_representative_hillslope(params,data_paths)
} else if (toolname == "define_run_printouts")  {
    define_run_printouts(params,data_paths)
}
else if (toolname == "define_run_printouts")  {
    define_run_printouts(params,data_paths)
}
else{
    # in any other case, the tool was invalid or not configured
    print(paste("[", Sys.time(), "] Either no TOOL_RUN environment variable available, or '", toolname, "' is not valid.\n", sep = "")) # nolint: line_length_linter.
}
