# ---- setup ----
dist_mat <- list()

coord_sites <- 
  data_sites %>% 
  column_to_rownames("site") %>% 
  dplyr::select(-number_samples)


# ---- geographic distance ----
# compute distance between locations (in meters)
dist_mat$geodist <- 
  geodist(coord_sites, measure = "geodesic") 

colnames(dist_mat$geodist) <- 
  rownames(dist_mat$geodist) <- 
  rownames(coord_sites)

dist_mat$geodist <-
  dist_mat$geodist %>% 
  as.dist()



# ---- distance by sea ----
library(marmap)
#load
Bathy <- 
  marmap::getNOAA.bathy(lon1 = -180, lon2 = 180,
                        lat1 = -30, lat2 = 30,
                        resolution = 10,
                        antimeridian = TRUE, # to be centered around pacific
                        keep = TRUE)

Bathy %>% saveRDS("intermediate/3_distance_metrics/bathymetry_10.RDS")


## transition matrix
trans <- 
  trans.mat(Bathy, min.depth=-5, max.depth=NULL)

trans %>% saveRDS("intermediate/3_distance_metrics/transition_matrix.RDS")


## distance by sea
dist_mat$seadist <- 
  lc.dist(trans, shift.lon(coord_sites), res = "dist")


# ---- export ----
dist_mat %>% 
  saveRDS("intermediate/3_distance_metrics/dist_geo.RDS")


