FROM rocker/rstudio:4.0.0
LABEL Maintainer="valavi.r@gmail.com"

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
    netcdf-bin \
    default-jre

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

# copy script and data
COPY R/ /home/rstudio/R
COPY data/ home/rstudio/data

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
RUN Rscript -e  "remotes::install_github('b0rxa/scmamp')"
RUN Rscript -e  "remotes::install_github('rvalavi/myspatial')"
RUN Rscript -e "remotes::install_version('gam', version = '1.20', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "remotes::install_version('gbm', version = '2.1.5', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('rJava')"

# copy maxent.jar file to dismo package
COPY maxent.jar /usr/local/lib/R/site-library/dismo/java/

# expose port
EXPOSE 8787
