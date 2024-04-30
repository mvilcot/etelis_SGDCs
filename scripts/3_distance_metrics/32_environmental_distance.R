
dist_mat <-
  readRDS("intermediate/3_distance_metrics/dist_geo.RDS")


# ---- load environmental predictors ----
## choose variables
variables <-
  c("BO2_chlomean_bdmean",
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



# ---- add buffer around stations ----
vect_stations <-
  data_stations %>% 
  dplyr::select(site, station, Longitude_approx, Latitude_approx) %>% 
  vect(geom = c("Longitude_approx", "Latitude_approx"), crs = "WGS84") %>%# convert to spatial vector 
  buffer(width = 50000) # add buffer


# ---- extract variables by station ----
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

envt_layers_Pcrop <- 
  envt_layers %>% 
  terra::rotate(left = FALSE) %>% # to put in 0/360 instead of -180/180
  terra::crop(extent(30, 220, -50, 50)) # crop to area of interest


## plot
pdf("results/3_distance_metrics/environmental_variables_bdmean_map2.pdf", width = 8, height = 5)
for(var in names(envt_layers_Pcrop)){
  gg <- 
    ggplot() +
    geom_spatraster(data = envt_layers_Pcrop[[var]]) +
    scale_fill_whitebox_c(palette = "muted") +
    geom_spatvector(data = vect_sites) +
    labs(fill = gsub("BO2_|_bdmean", "", var)) +
    theme_light() +
    theme(panel.grid.major = element_blank(),
          legend.position="bottom") 
  print(gg)
}
dev.off()





# ---- PCA ----
envt_var <- 
  envt_site %>% 
  column_to_rownames("site")

library(corrplot)
corrplot(cor(envt_var))

# run pca
pca_env <- dudi.pca(df = envt_var, scannf = F, nf = 5)
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

#### {FIGURE S8} ####
# plot
LABELS <- gsub("_", " ", rownames(pca_env$li))
LABELS <- gsub('W Australia', 'Western Australia', LABELS)
gg <-
  ggplot(pca_env$li, aes(x=Axis1, y=Axis2, color=rownames(pca_env$li))) +
  geom_point() +
  labs(x=labels[1], y=labels[2]) +
  ggrepel::geom_text_repel(label = LABELS) +
  scale_color_manual(values = color_perso) +
  theme_light() +
  theme(legend.position="none")
gg

ggsave("results/3_distance_metrics/_S8_PCA_environmental_variables_bdmean_absolute.png",
       gg,
       height = 5.5, width = 6, units = "in", dpi = 300)



# ---- environmental distance ----
# OR more directly, euclidian distance between scaled values
dist_mat$environment <-
  vegdist(envt_var, method = "euclidean", na.rm = TRUE) 


# export
dist_mat %>% 
  saveRDS("intermediate/3_distance_metrics/dist_geo_envtbdmean.RDS")

