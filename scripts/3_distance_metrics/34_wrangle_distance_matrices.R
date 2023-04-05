
dist_mat <-
  readRDS("intermediate/3_distance_metrics/dist_geo_envt_res.RDS")



## ---- pivot longer distance matrix ----
# melt_dist <- 
#   dist_mat %>% 
#   lapply(function(i, n){melt.dist(dist = i, metric = names(i))})

dist_melt <- list()
for (i in 1:length(dist_mat)){
  dist_melt[[i]] <- melt.dist(distmat = dist_mat[[i]], metric = names(dist_mat)[[i]])
}


## ---- merge into a table ----
dist_merge <- 
  dist_melt %>%
  reduce(full_join, by = c("site1", "site2"))  

pairs(dist_merge[3:6])


## ---- save ----
dist_merge %>% 
  write.csv("results/3_distance_metrics/dist_geo_envt_res.csv", 
            row.names = F, quote = F)


