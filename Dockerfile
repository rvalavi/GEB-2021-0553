FROM rocker/rstudio:4.0.0
LABEL Maintainer="valav.r@gmail.com"

# install goespatial libraries
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    gcc \
    libgdal-dev \
    libgeos-dev \
    libgeotiff-dev \
    libproj-dev \
    libudunits2-dev \
    libhdf5-dev \
    libnetcdf-dev \
    netcdf-bin

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

# change working directory
WORKDIR /project
# recover R packages
COPY renv.lock renv.lock
RUN Rscript -e "install.packages('renv')"
RUN Rscript -e "renv::consent(provided = TRUE)"
RUN Rscript -e "renv::restore()"
# install GitHub packages
RUN Rscript -e "install.packages('remotes')"
RUN Rscript -e  "remotes::install_github('meeliskull/prg/R_package/prg')"
RUN Rscript -e  "remotes::install_github('rvalavi/myspatial')"

# expose port
EXPOSE 8787
