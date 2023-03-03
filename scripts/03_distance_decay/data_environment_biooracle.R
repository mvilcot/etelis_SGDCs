# script to extract environmental variables from the bio-oracle database

# set ----
# packages
library(tidyverse)
library(readxl)
library(terra)
library(sdmpredictors)
level = "site"


layers_velocity <- 
  list_layers() %>%
  filter(dataset_code == "Bio-ORACLE") %>% 
  filter(grepl("Current velocity", name))
write.csv(layers_velocity, "intermediate/00_sampling_sites/layers_bio-oracle_velocity.csv",
          quote = F, row.names = F)


# variables
variables <-
  c("BO_chlomean",
    "BO_dissox",
    "BO_nitrate",
    "BO_ph",
    "BO_salinity",
    "BO_sstmean",
    "BO_sstmax",
    "BO_sstmin",
    "BO21_tempmean_bdmean",
    "BO22_curvelltmax_bdmean")

# load ----
# environmental predictors
data_variables <-
  sdmpredictors::load_layers(variables)



# wrangle data ----
# convert rasters to spatRaster
data_variables <- rast(data_variables)

# get mean coordinates by location
data_coord <- 
  data_sites %>% 
  group_by(.data[[level]]) %>% 
  summarise(longitude=mean(Longitude_approx),
            latitude=mean(Latitude_approx)) %>% 
  column_to_rownames(level)

vect_coord <-
  data_coord %>% 
  vect(geom = c("longitude", "latitude"))     # convert to spatial vector 

# extract variables by site ----
site_variables <- 
  extract(data_variables,
          vect_coord,
          xy = TRUE,
          ID = FALSE,
          bind = TRUE)

site_variables <- 
  as.data.frame(site_variables) %>% 
  rename(longitude = "x", latitude = "y")
  cbind(data_coord)

# export ----
write_csv(site_variables,
          "./data/environment/bio_oracle/bio_oracle_variables.csv")





