comm_delin <- "taxonomy"
dist_merge <- 
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_merge$site1 <- factor(dist_merge$site1, levels = levels(data_sites$site))
dist_merge$site2 <- factor(dist_merge$site2, levels = levels(data_sites$site))
dist_merge

dist_mergeSUB <- 
  dist_merge %>% 
  mat.subset("Seychelles")


# ---- beta GD heatmap ----
metricGD <- "D.Jost"

gg1 <- 
  ggplot(data = dist_merge, aes(site2, site1, fill = .data[[metricGD]]))+ 
  geom_tile()+ 
  scale_fill_viridis(direction = -1, option = "rocket", begin = 0.05, end = 0.95) +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) + 
  coord_fixed()


gg2 <- 
  ggplot(data = dist_mergeSUB, aes(site2, site1, fill = .data[[metricGD]]))+ 
  geom_tile()+ 
  scale_fill_viridis(direction = -1, option = "rocket", begin = 0.05, end = 0.95) +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) + 
  coord_fixed()


gg1 + gg2 +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  patchwork::plot_annotation(title = "Etelis coruscans")

ggsave(paste0("results/3_distance_metrics/heatmap_GD_", metricGD, ".png"),
       width = 10, height = 4)




# ---- beta SD heatmap ----
comm <- "Lutjanidae"

metricSDcomm <- paste0(comm, ".beta.jac")
gg1 <- 
  ggplot(data = dist_merge, aes(site2, site1, fill = .data[[metricSDcomm]]))+ 
  geom_tile()+ 
  scale_fill_viridis(direction = -1, option = "mako", begin = 0.05, end = 0.95) +
  labs(fill = "beta.jac") +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) + 
  coord_fixed()


metricSDcomm <- paste0(comm, ".beta.jtu")
gg2 <- 
  ggplot(data = dist_merge, aes(site2, site1, fill = .data[[metricSDcomm]]))+ 
  geom_tile()+ 
  scale_fill_viridis(direction = -1, option = "mako", begin = 0.05, end = 0.95) +
  labs(fill = "beta.jtu") +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) + 
  coord_fixed()


gg1 + gg2 +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  patchwork::plot_annotation(title = comm)

ggsave(paste0("results/3_distance_metrics/heatmap_SD_", comm, ".png"),
       width = 10, height = 4)



# ---- seadist heatmap ----

gg <-  
  ggplot(data = dist_merge, aes(site2, site1, fill = seadist))+ 
  geom_tile()+
  scale_fill_gradient(low = "#DEEBF7", high = "#3182BD") +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) + 
  coord_fixed()


gg

ggsave(paste0("results/3_distance_metrics/heatmap_seadist.png"),
       width = 5, height = 4)




# ---- environment heatmap ----

brewer.pal(n = 3, name = "Greens")

gg <-  
  ggplot(data = dist_merge, aes(site2, site1, fill = environment))+ 
  geom_tile()+ 
  scale_fill_gradient(low = "#E5F5E0", high = "#31A354") +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) + 
  coord_fixed()


gg

ggsave(paste0("results/3_distance_metrics/heatmap_environment.png"),
       width = 5, height = 4)

