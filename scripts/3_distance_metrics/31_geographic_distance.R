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
Bathy <- readRDS("intermediate/0_sampling_design/_archive/bathymetry.RDS")

# compute distance by sea
trans <- 
  trans.mat(Bathy, min.depth=-5, max.depth=NULL)

#### >>>>> !!! error when too high resolution for Guam ############
dist_mat$seadist <- 
  lc.dist(trans, shift.lon(coord_sites), res = "dist")

get.depth(Bathy, shift.lon(coord_sites), locator = F)

# ---- export ----
dist_mat %>% 
  saveRDS("intermediate/3_distance_metrics/dist_geo2.RDS")
# dist_mat2 <-
#   readRDS("intermediate/3_distance_metrics/dist_geo.RDS")


