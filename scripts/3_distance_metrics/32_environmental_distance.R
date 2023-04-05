
dist_mat <-
  readRDS("intermediate/3_distance_metrics/dist_geo.RDS")


# ## ---- current velocity ----
# layers_velocity <- 
#   list_layers() %>%
#   filter(dataset_code == "Bio-ORACLE") %>% 
#   filter(grepl("Current velocity", name))
# write.csv(layers_velocity, "intermediate/0_sampling_design/layers_bio-oracle_velocity.csv",
#           quote = F, row.names = F)


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

envt_layers <-
  sdmpredictors::load_layers(variables)


## ---- wrangle data ----
# convert rasters to spatRaster
envt_layers <- 
  rast(envt_layers)

# convert to spatial vector 
vect_coord <-
  coord_site %>% 
  vect(geom = c("longitude", "latitude"), crs = "WGS84")

## ---- extract variables by site ----
envt_site <- 
  terra::extract(envt_layers,
                 vect_coord,
                 xy = TRUE,
                 ID = FALSE,
                 bind = TRUE)

envt_site <- 
  as.data.frame(envt_site) %>% 
  cbind(coord_site) 



## ---- export ----
envt_site %>% 
  write.csv("intermediate/3_distance_metrics/bio_oracle_variables.csv",
            row.names = T, quote = F)



# ## ---- plot ----
# # Generate a nice color ramp and plot the map
# my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
# plot(envt_layers$BO_salinity, col=my.colors(1000))
# title(cex.sub = 1.25, sub = "Maximum temperature at the sea bottom (ºC)")
# 




## ---- Matrix environmental variables ----
# read variables 

# keep environmental variables
envt_var <- 
  envt_site %>% 
  dplyr::select(-c(x, y, longitude, latitude))

# run pca
pca_env <- dudi.pca(df = envt_var, scannf = FALSE, nf = 5) 

# compute euclidian distance between coordinates
dist_mat$environment <- 
  vegdist(pca_env$li, method = "euclidean", na.rm = TRUE) # euclidean distance based on the 3 PCA axes

# extract pca coordinates
pca_coord <- 
  pca_env$li %>% 
  rownames_to_column("site")

# # analyse pca
# library("factoextra")
# pca_eig <- get_eigenvalue(pca_env)
# fviz_eig(pca_env, addlabels = TRUE, ylim = c(0, 50))
# fviz_pca_var(pca_env, col.var = "black")



## ---- export ----
dist_mat %>% 
  saveRDS("intermediate/3_distance_metrics/dist_geo_envt.RDS")



