# script to extract environmental variables from the bio-oracle database

## ---- parameters ----
level = "site"


layers_velocity <- 
  list_layers() %>%
  filter(dataset_code == "Bio-ORACLE") %>% 
  filter(grepl("Current velocity", name))
write.csv(layers_velocity, "intermediate/0_sampling_design/layers_bio-oracle_velocity.csv",
          quote = F, row.names = F)


## ---- load environmental predictors ----
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

data_variables <-
  sdmpredictors::load_layers(variables)



## ---- wrangle data ----
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

## ---- extract variables by site ----
site_variables <- 
  extract(data_variables,
          vect_coord,
          xy = TRUE,
          ID = FALSE,
          bind = TRUE)

site_variables <- 
  as.data.frame(site_variables) %>% 
  # dplyr::rename(longitude = "x", latitude = "y") %>% 
  cbind(data_coord) %>% 
  rownames_to_column("site")


## ---- export ----
write.csv(site_variables, 
          "intermediate/3_distance_decay/bio_oracle_variables.csv", 
          row.names = F, quote = F)


## ---- plot ----
# Generate a nice color ramp and plot the map
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
plot(data_variables$BO_salinity, col=my.colors(1000))
title(cex.sub = 1.25, sub = "Maximum temperature at the sea bottom (ºC)")

