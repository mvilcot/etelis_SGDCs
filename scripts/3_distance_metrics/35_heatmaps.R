comm_delin <- "taxonomy"
dist_merge <- 
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_merge$site1 <- factor(dist_merge$site1, levels = levels(data_sites$site))
dist_merge$site2 <- factor(dist_merge$site2, levels = levels(data_sites$site))
dist_merge

dist_mergeSUB <- 
  dist_merge %>% 
  mat.subset("Seychelles")

names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")



# great labels ----
dist_merge <- 
  dist_merge %>% 
  mutate(site1 = gsub('_', ' ', site1)) %>% 
  mutate(site2 = gsub('_', ' ', site2)) %>% 
  mutate(site1 = gsub('W Australia', 'Western Australia', site1)) %>% 
  mutate(site2 = gsub('W Australia', 'Western Australia', site2)) %>% 
  mutate(site1 = factor(site1, levels = LABELS)) %>% 
  mutate(site2 = factor(site2, levels = LABELS))
  
dist_mergeSUB <- 
  dist_mergeSUB %>% 
  mutate(site1 = gsub('_', ' ', site1)) %>% 
  mutate(site2 = gsub('_', ' ', site2)) %>% 
  mutate(site1 = gsub('W Australia', 'Western Australia', site1)) %>% 
  mutate(site2 = gsub('W Australia', 'Western Australia', site2)) %>% 
  mutate(site1 = factor(site1, levels = LABELS[-1])) %>% 
  mutate(site2 = factor(site2, levels = LABELS[-1]))

gglist <- list()


# ---- seadist heatmap ----
gg1 <-  
  ggplot(data = dist_merge, aes(site2, site1, fill = seadist))+ 
  geom_tile()+
  scale_fill_gradient(low = "#DEEBF7", high = "#3182BD") +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) + 
  ggtitle("Distance by sea") +
  labs(fill = "km") +
  coord_fixed()



# ---- environment heatmap ----
gg2 <-  
  ggplot(data = dist_merge, aes(site2, site1, fill = environment))+ 
  geom_tile()+ 
  scale_fill_gradient(low = "#E5F5E0", high = "#31A354") +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) +
  ggtitle("Environmental distance") +
  labs(fill = "") +
  coord_fixed()

gglist[["envt"]] <- 
  gg1 + gg2



# ---- beta GD heatmap ----
metricGD <- "Fst"

gg1 <-  
  ggplot(data = dist_merge, aes(site2, site1, fill = .data[[metricGD]]))+ 
  geom_tile()+ 
  scale_fill_viridis(direction = -1, option = "rocket", begin = 0.05, end = 0.95) +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) + 
  ggtitle("Etelis coruscans")+
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
  ggtitle("Etelis coruscans") +
  coord_fixed()

gglist[[metricGD]] <-
  gg1 + gg2 



# ---- beta SD heatmap ----
for (comm in names_communities){
  
  metricSDcomm <- paste0(comm, ".beta.jac")
  gg1 <- 
    ggplot(data = dist_merge, aes(site2, site1, fill = .data[[metricSDcomm]]))+ 
    geom_tile()+ 
    scale_fill_viridis(direction = -1, option = "mako", begin = 0.05, end = 0.95) +
    labs(fill = "βjac") +
    labs(x = "", y = "") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          panel.grid = element_blank(), 
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")) + 
    ggtitle(comm) +
    coord_fixed()
  
  
  metricSDcomm <- paste0(comm, ".beta.jtu")
  gg2 <- 
    ggplot(data = dist_merge, aes(site2, site1, fill = .data[[metricSDcomm]]))+ 
    geom_tile()+ 
    scale_fill_viridis(direction = -1, option = "mako", begin = 0.05, end = 0.95) +
    labs(fill = "βjtu") +
    labs(x = "", y = "") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          panel.grid = element_blank(), 
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")) + 
    ggtitle(comm) +
    coord_fixed()
  
  
  gglist[[comm]] <-
    gg1 + gg2 +
    patchwork::plot_annotation(title = comm)
  
}


#### {FIGURE S4} ####
# all heatmaps together ----
ggall <- 
  wrap_plots(gglist, ncol = 1) + 
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
  
ggsave(paste0("results/3_distance_metrics/_S4_heatmaps_all.png"),
       plot = ggall,
       width = 8, height = 20, dpi = 500)


