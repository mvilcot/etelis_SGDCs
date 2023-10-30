# ---- read distance matrix ----
level = "site"
comm_delin = "taxonomy"
metric <- "Fst" # Lutjanidat.beta.jtu

dist_merge <- read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))
dist_mat <- readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))

explanatory <- 
  read_csv("intermediate/3_distance_metrics/bio_oracle_variables_bdmean_stationbuffer_siteaverage.csv") %>% 
  left_join(data_sites, by = "site")


# # ---- subset Seychelles ----
# explanatory <- mat.subset(explanatory, "Seychelles")
# dist_mat[[metric]] <- mat.subset(dist_mat[[metric]], "Seychelles")


# ---- Run RDA ----

dbrda_run  <- 
  capscale(
    dist_mat[[metric]] ~ 
      # longitude + latitude,
      BO2_chlomean_bdmean +
      BO2_dissoxmean_bdmean +
      BO2_nitratemean_bdmean +
      BO2_tempmean_bdmean +
      BO2_salinitymean_bdmean +
      Condition(longitude + latitude),
    data = explanatory)

# Get the inertia statistics                   
dbrda_run

# Test the significance of the model
anova.cca(dbrda_run, permutations = 9999)

# Calculate the adjusted-R²
RsquareAdj(dbrda_run)


# ---- Plot ----
# automatic
# plot(dbrda_run)

# get data
dbrda_sites <- 
  as.data.frame(vegan::scores(dbrda_run)$sites) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("site") %>% 
  left_join(data_sites, by = "site")

dbrda_var <- 
  as.data.frame(dbrda_run[["CCA"]][["biplot"]][, 1:2]) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("variable")
  
ggplot() +
  # geom_encircle(aes(group = Protection,linetype = Protection,fill= Protection), s_shape = 1, expand = 0,
  #               alpha = 0.4, show.legend = FALSE) + # hull area 
  geom_point(data = dbrda_sites, aes(x= CAP1, y = CAP2, color = site)) +
  geom_text_repel(data = dbrda_sites, aes(x= CAP1, y = CAP2, color = site),
                  label = dbrda_sites$site, bg.color = "grey70", bg.r = 0.02) +
  scale_color_manual(values = color_perso) +
  geom_segment(data= dbrda_var, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc")), show.legend = F) +
  geom_label_repel(data= dbrda_var, 
                   aes(x= CAP1, y=CAP2, fontface=3),
                   label = dbrda_var$variable,
                   label.size = NA,
                   size = 4,
                   fill = NA) +
  ggtitle(metric) +
  theme_light() +
  theme(legend.position="none")

