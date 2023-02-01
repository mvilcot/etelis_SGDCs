# script to extract environmental variables from the bio-oracle database

# set ----
# packages
library(tidyverse)
library(readxl)
library(terra)
library(sdmpredictors)

# variables
variables <-
  c("BO_chlomean",
    "BO_dissox",
    "BO_nitrate",
    "BO_ph",
    "BO_salinity",
    "BO_sstmean",
    "BO_sstmax",
    "BO_sstmin")

# load ----
# environmental predictors
data_variables <-
  load_layers(variables)

# sample meta data
data_sample <-
  read_excel("./data/sample_meta/data_sample.xlsx",
             guess_max = 10000)

data_rad <-
  read_excel("./data/sample_meta/data_RAD.xlsx")

# wrangle data ----
# convert rasters to spatRaster
data_variables <- rast(data_variables)

# sample coordinates
data_coords <-
  left_join(data_rad,data_sample) %>%       # filter only rad samples
  dplyr::select(c(longitude,latitude,sample_site)) %>%  # retain only coord data
  distinct() %>%                            # remove duplicates
  vect(geom=c("longitude", "latitude"))     # convert to spatial vector 

# extract variables by site ----
site_variables <- 
  extract(data_variables,
          data_coords,
          xy = TRUE,
          ID = FALSE,
          bind = TRUE)

site_variables <- 
  as.data.frame(site_variables)

# export ----
write_csv(site_variables,
          "./data/environment/bio_oracle/bio_oracle_variables.csv")





