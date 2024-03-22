
dist_mat <-
  readRDS("intermediate/3_distance_metrics/dist_geo.RDS")


# ---- load environmental predictors ----
## choose variables
variables <-
  c(# "BO2_ph", # does not exist for bottom depth
    "BO2_chlomean_bdmean",
    "BO2_dissoxmean_bdmean",
    "BO2_nitratemean_bdmean",
    "BO2_salinitymean_bdmean",
    "BO2_tempmean_bdmean",
    "BO2_tempmax_bdmean",
    "BO2_tempmin_bdmean"
  )

## load online
envt_layers <-
  sdmpredictors::load_layers(variables,
                             datadir = "intermediate/3_distance_metrics/sdmpredictors")

envt_layers <-
  rast(envt_layers)

# ## OR load from downloaded file
# envt_layers <- 
#   rast(list.files("intermediate/3_distance_metrics/bio-oracle", full.names = T))



# ---- add buffer around stations ----
vect_stations <-
  data_stations %>% 
  dplyr::select(site, station, Longitude_approx, Latitude_approx) %>% 
  vect(geom = c("Longitude_approx", "Latitude_approx"), crs = "WGS84") %>%# convert to spatial vector 
  buffer(width = 50000) # add buffer


# ---- extract variables by station ----
## >>>> ERROR MHI AND NEW CALEDONIA if too small buffer... ----
envt_station <- 
  terra::extract(envt_layers,
                 vect_stations,
                 xy = TRUE,
                 ID = FALSE,
                 bind = TRUE,
                 fun = mean, # average across station buffer
                 na.rm = TRUE) %>% 
  as_tibble()


## export
envt_station %>% 
  write_csv("intermediate/3_distance_metrics/bio_oracle_variables_bdmean_stationbuffer.csv")


# ---- average by site ----
## >>>>> TRY TO AVERAGE ACROSS ALL STATION BUFFER DIRECTLY ##########
envt_site <-
  envt_station %>%
  pivot_longer(-c(site, station), names_to = "variable") %>% 
  group_by(site, variable) %>% 
  summarise(value = mean(value), .groups = "keep") %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  ungroup()


## export
envt_site %>% 
  write_csv("intermediate/3_distance_metrics/bio_oracle_variables_bdmean_stationbuffer_siteaverage.csv")


# ---- scale data ----
envt_site_scaled <-
  envt_site %>%
  # mutate(across(where(is.numeric), log)) %>%
  column_to_rownames("site") %>%
  scale() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("site") ## the rownames does not appear but they are still stored...

# pivot longer
envt_site_long <-
  envt_site %>% 
  pivot_longer(cols = -site, names_to = "variable", values_to = "value")

envt_site_long_scaled <-
  envt_site_scaled %>%
  pivot_longer(cols = -site, names_to = "variable", values_to = "value")

# relevel
envt_site_long$site <- 
  envt_site_long$site %>%  
  ordered(levels = unique(data_stations[order(data_stations$order),][["site"]]))

envt_site_long_scaled$site <-
  envt_site_long_scaled$site %>%
  ordered(levels = unique(data_stations[order(data_stations$order),][["site"]]))



# ---- boxplot ----
ggplot(envt_site_long, aes(site, value, fill = site, color = site)) +
  geom_boxplot() +
  facet_wrap( ~ variable, ncol=1, scales = "free_y") +
  scale_color_manual(values = color_perso) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(width = 5, height = 10, 
       filename = paste0("results/3_distance_metrics/environmental_variables_bdmean_absolute.png"))

ggplot(envt_site_long_scaled, aes(site, value, fill = site, color = site)) +
  geom_boxplot() +
  facet_wrap( ~ variable, ncol=1, scales = "free_y") +
  scale_color_manual(values = color_perso) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(width = 5, height = 10,
       filename = paste0("results/3_distance_metrics/environmental_variables_bdmean_scaled.png"))



# ---- map ----
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

## center to pacific
vect_sites <-
  data_sites %>% 
  vect(geom = c("longitude", "latitude"), crs = "WGS84") %>% 
  terra::rotate(left = FALSE) # to put in 0/360 instead of -180/180

# vect_sites <- # do the same thing, manually
#   vect(shift.lon(data_sites), geom=c("longitude", "latitude"), 
#        crs = "EPSG:4326")

envt_layers_Pcrop <- 
  envt_layers %>% 
  terra::rotate(left = FALSE) %>% # to put in 0/360 instead of -180/180
  terra::crop(extent(30, 220, -50, 50)) # crop to area of interest


## plot
pdf("results/3_distance_metrics/environmental_variables_bdmean_map2.pdf", width = 8, height = 5)
for(var in names(envt_layers_Pcrop)){
  # plot(envt_layers_Pcrop[[var]], col=my.colors(1000))
  # plot(vect_sites, add = T, col = "grey20")
  # title(cex.sub = 1.25, main = var)
  
  gg <- 
    ggplot() +
    geom_spatraster(data = envt_layers_Pcrop[[var]]) +
    scale_fill_whitebox_c(palette = "muted") +
    geom_spatvector(data = vect_sites) +
    labs(fill = gsub("BO2_|_bdmean", "", var)) +
    # scale_fill_viridis_c() +
    # scale_fill_distiller(palette = "YlGnBu") +
    theme_light() +
    theme(panel.grid.major = element_blank(),
          legend.position="bottom") 
  print(gg)
}
dev.off()





# ---- PCA ----
# read variables 

# keep environmental variables
# envt_var <- 
#   envt_site %>% 
#   dplyr::select(-c(x, y, longitude, latitude))

envt_var <- 
  envt_site_scaled %>% 
  column_to_rownames("site")

library(corrplot)
corrplot(cor(envt_var))

# run pca
pca_env <- dudi.pca(df = envt_var, scannf = FALSE, nf = 5)
pca_env$li %>%
  rownames_to_column("site") %>%
  write_csv("intermediate/3_distance_metrics/PCA_environment_axis_site.csv")

# analyse pca
library(factoextra)
pca_eig <- get_eigenvalue(pca_env)
fviz_eig(pca_env, addlabels = TRUE)
fviz_pca_var(pca_env)
fviz_contrib(pca_env, choice = "var", axes = 1, top = 10)

# eigenvalues
percent_explained <- pca_env$eig / sum(pca_env$eig) * 100

# set axis labels
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("PC1 ({pretty_pe[1]}%)"),
            glue("PC2 ({pretty_pe[2]}%)"))

# plot
gg <-
  ggplot(pca_env$li, aes(x=Axis1, y=Axis2, color=rownames(pca_env$li))) +
  geom_point() +
  labs(x=labels[1], y=labels[2]) +
  ggrepel::geom_text_repel(label = rownames(pca_env$li)) +
  scale_color_manual(values = color_perso) +
  theme_light() +
  theme(legend.position="none")
ggsave("results/3_distance_metrics/PCA_environmental_variables_bdmean_scaled.png",
       gg,
       height = 7, width = 7, units = "in", dpi = 300)


# # eigenvalues
# percent_explained <- pca_env$eig / sum(pca_env$eig) * 100
# 


# ---- environmental distance ----
# # euclidian distance from the first 3 PCA axes
# dist_mat$environment <-
#   vegdist(pca_env$li, method = "euclidean", na.rm = TRUE) #

# OR more directly, euclidian distance between scaled values
dist_mat$environment <-
  vegdist(envt_var, method = "euclidean", na.rm = TRUE) 


# export
dist_mat %>% 
  saveRDS("intermediate/3_distance_metrics/dist_geo_envtbdmean.RDS")




# ---- *** DRAFTS ----
## *** current velocity ----
# layers_velocity <- 
#   list_layers() %>%
#   filter(dataset_code == "Bio-ORACLE") %>% 
#   filter(grepl("Current velocity", name))
# write.csv(layers_velocity, "intermediate/0_sampling_design/layers_bio-oracle_velocity.csv",
#           quote = F, row.names = F)

