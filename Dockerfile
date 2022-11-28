# Pull any base image that includes R
FROM r-base:4.2.0

# the parameter parsing function always needs the rjson and yaml packages
RUN R -e "install.packages(c('jsonlite', 'yaml'))"

# install json2aRgs v0.3.0 to parse parameters from /in/parameters.json
RUN wget -q https://cran.r-project.org/src/contrib/json2aRgs_0.3.0.tar.gz
RUN R CMD INSTALL json2aRgs_0.3.0.tar.gz

# remove json2aRgs tarball
RUN rm json2aRgs_0.3.0.tar.gz

# install Catflow-R-Package dependencies
RUN R -e "install.packages(c('deSolve', 'RColorBrewer', 'zoo', 'xts'))"

# create the tool input structure
RUN mkdir /in
COPY ./in /in
RUN mkdir /out
RUN mkdir /src
COPY ./src /src

# download latest version of Catflow-R-Package from github, untar and rename
RUN wget -qO- https://github.com/CATFLOW/Catflow-R-Package/archive/refs/heads/main.tar.gz | tar xz -C src/ && mv /src/Catflow-R-Package-main /src/Catflow-R-Package

# install Catflow-R-Package from source
RUN R -e "install.packages('/src/Catflow-R-Package/Catflow', repos = NULL, type = 'source')"

# delete Catflow-R-Pacakge source code
RUN rm -rf /src/Catflow-R-Package

WORKDIR /src
CMD ["Rscript", "run.R"]
