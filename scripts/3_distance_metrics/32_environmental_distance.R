
dist_mat <-
  readRDS("intermediate/3_distance_metrics/dist_geo.RDS")


# ---- load environmental predictors ----
variables <-
  c("BO2_chlomean_bdmean",
    "BO2_dissoxmean_bdmean",
    "BO2_nitratemean_bdmean",
    # "BO2_ph", # does not exist for bottom depth
    "BO2_salinitymean_bdmean",
    "BO2_tempmean_bdmean",
    "BO2_tempmax_bdmean",
    "BO2_tempmin_bdmean"
  )

envt_layers <-
  sdmpredictors::load_layers(variables)


# ---- setup data ----
# convert rasters to spatRaster
envt_layers <- 
  rast(envt_layers)

# convert to spatial vector 
vect_coord <-
  coord_site %>% 
  vect(geom = c("longitude", "latitude"), crs = "WGS84")


# ---- extract variables by site ----
envt_site <- 
  terra::extract(envt_layers,
                 vect_coord,
                 xy = TRUE,
                 ID = FALSE,
                 bind = TRUE)

envt_site <- 
  as.data.frame(envt_site) %>% 
  cbind(coord_site) 


## export
envt_site %>% 
  write.csv("intermediate/3_distance_metrics/bio_oracle_variables_bdmean.csv",
            row.names = T, quote = F)


# ---- scale variables ----
envt_site_scaled <- 
  envt_site %>%
  dplyr::select(all_of(variables)) %>% 
  log() %>% ################# CHECK IF LOG NEEDED (cf DiBattista et al. 2020)
  scale() %>% 
  as.data.frame() %>% 
  rownames_to_column("site")

# pivot longer
envt_site2 <-
  envt_site %>% 
  dplyr::select(all_of(variables)) %>% 
  rownames_to_column("site") %>% 
  pivot_longer(cols = variables, names_to = "variable", values_to = "value")

envt_site2_scaled <-
  envt_site_scaled %>% 
  pivot_longer(cols = variables, names_to = "variable", values_to = "value")

# relevel
envt_site2$site <- 
  envt_site2$site %>%  
  ordered(levels = unique(data_sites[order(data_sites$order),][["site"]]))

envt_site2_scaled$site <- 
  envt_site2_scaled$site %>%  
  ordered(levels = unique(data_sites[order(data_sites$order),][["site"]]))



# ---- boxplot ----
ggplot(envt_site2, aes(site, value, fill = site, color = site)) +
  geom_boxplot() +
  facet_wrap( ~ variable, ncol=1, scales = "free_y") +
  scale_color_manual(values = color_perso) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(width = 5, height = 10, 
       filename = paste0("results/3_distance_metrics/environmental_variables_absolute.png"))

ggplot(envt_site2_scaled, aes(site, value, fill = site, color = site)) +
  geom_boxplot() +
  facet_wrap( ~ variable, ncol=1, scales = "free_y") +
  scale_color_manual(values = color_perso) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(width = 5, height = 10, 
       filename = paste0("results/3_distance_metrics/environmental_variables_scaled.png"))



# ---- map plot----
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

coord_vect <- vect(shift.lon(coord_site), geom=c("longitude", "latitude"), crs = "EPSG:4326")

envt_layers_360 <- 
  terra::rotate(envt_layers, left = FALSE)

envt_layers_crop <- 
  terra::crop(envt_layers_360, extent(30, 220, -35, 35))

pdf("results/3_distance_metrics/environmental_variables_plot.pdf", width = 12, height = 5)
for(var in variables){
  plot(envt_layers_crop[[var]], col=my.colors(1000))
  plot(coord_vect, add = T, col = "grey20")
  title(cex.sub = 1.25, main = var)
}
dev.off()


# library(tidyterra)
# ggplot() +
#   geom_spatraster(data = envt_layers_crop[[variables[1]]]) +
#   scale_fill_whitebox_c(palette = "muted") +
#   theme_light()




# ---- PCA ----
# read variables 

# keep environmental variables
# envt_var <- 
#   envt_site %>% 
#   dplyr::select(-c(x, y, longitude, latitude))

envt_var <- 
  envt_site_scaled %>% 
  column_to_rownames("site")

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

ggplot(pca_coord, aes(x=Axis1, y=Axis2, color = site)) +
  geom_point() +
  # labs(x=labels[1], y=labels[2]) +
  geom_text(aes(label = site)) +
  scale_color_manual(values = color_perso) +
  theme(legend.position="none")
ggsave("results/3_distance_metrics/environmental_variables_scaled_pca.png", 
       height = 7, width = 7, units = "in", dpi = 300)


## export
dist_mat %>% 
  saveRDS("intermediate/3_distance_metrics/dist_geo_envt.RDS")




# ---- *** DRAFTs ----
## *** current velocity ----
# layers_velocity <- 
#   list_layers() %>%
#   filter(dataset_code == "Bio-ORACLE") %>% 
#   filter(grepl("Current velocity", name))
# write.csv(layers_velocity, "intermediate/0_sampling_design/layers_bio-oracle_velocity.csv",
#           quote = F, row.names = F)

