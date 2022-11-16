# load the parameters parser
source("get_parameters.R")

# get the call parameters for the tool
params <- get_parameters()

# check which tool to run
toolname <- tolower(Sys.getenv("TOOL_RUN"))

if (toolname == "") {
    toolname <- "foobar"
}

# Switch for the different tools available in this package
if (toolname == "foobar") {
    # RUN the tool here and create the output in /out
    
    # prints to STDOUT.log
    print("You have tried to run the tool 'foobar'.
This tool is the template tool without any functionality.
Please implement another tool or select the tool you 
have already implemented.")
    print(params)
} else {
    # in any other case, the tool was invalid or not configured
    print(paste("[", Sys.time(), "] Either no TOOL_RUN environment variable available, or '", toolname, "' is not valid.\n", sep = ""))
}
