# Pull any base image that includes R
FROM rocker/geospatial:4.4.1

# the parameter parsing function always needs the rjson and yaml packages
RUN R -e "install.packages(c('jsonlite', 'yaml'))"

# install json2aRgs latest version from github
RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_github("VForWaTer/json2aRgs")'

# install Catflow-R-Package dependencies
RUN R -e "install.packages(c('deSolve', 'RColorBrewer', 'zoo', 'xts'))"

# create the tool input structure
RUN mkdir /in
COPY ./in /in
RUN mkdir /out
RUN mkdir /src
COPY ./src /src

RUN R -e "install.packages(c('raster', 'WVPlots', 'sigr', 'gridExtra', 'mgcv', 'ggplot2', 'wrapr', 'dplyr', 'ggbeeswarm', 'ggstatsplot'))"

# install Catflow-R-Package from source
RUN R -e "install.packages('/src/Catflow_0.998_v2.tar.gz', repos = NULL, type = 'source')"

# install representative hillslope dependencies
# Download rgdal package archive, untar, and rename
RUN wget -qO- https://cran.r-project.org/src/contrib/Archive/rgdal/rgdal_1.6-7.tar.gz | tar xz -C src/ 

# Install rgdal package from source
RUN R -e "install.packages('/src/rgdal', repos = NULL, type = 'source')" #


WORKDIR /src
CMD ["Rscript", "run.R"]
