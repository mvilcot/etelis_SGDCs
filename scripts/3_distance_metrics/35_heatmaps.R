comm_delin <- "taxonomy"
dist_merge <- 
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_merge$site1 <- factor(dist_merge$site1, levels = levels(data_sites$site))
dist_merge$site2 <- factor(dist_merge$site2, levels = levels(data_sites$site))
dist_merge

dist_mergeSUB <- 
  dist_merge %>% 
  mat.subset("Seychelles")

names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria/misc", "Teleostei")


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



# ggsave(paste0("results/3_distance_metrics/heatmap_seadist.png"),
#        width = 5, height = 4)




# ---- environment heatmap ----

# RColorBrewer::brewer.pal(n = 3, name = "Greens")

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


# ggsave(paste0("results/3_distance_metrics/heatmap_environment.png"),
#        width = 5, height = 4)

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
  ggtitle("Etelis coruscans")+
  coord_fixed()

gglist[[metricGD]] <-
  gg1 + gg2 

# gg1 + gg2 +
#   patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
#   patchwork::plot_annotation(title = "Etelis coruscans")
# 
# ggsave(paste0("results/3_distance_metrics/heatmap_GD_", metricGD, ".png"),
#        width = 10, height = 4)




# ---- beta SD heatmap ----

for (comm in names_communities){
  
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
    ggtitle(comm) +
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
    ggtitle(comm) +
    coord_fixed()
  
  
  gglist[[comm]] <-
    gg1 + gg2 +
    patchwork::plot_annotation(title = comm)
  
}


# gg1 + gg2 +
#   patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
#   patchwork::plot_annotation(title = comm)
# 
# ggsave(paste0("results/3_distance_metrics/heatmap_SD_", comm, ".png"),
#        width = 10, height = 4)



# all maps together ----

ggall <- 
  wrap_plots(gglist, ncol = 1) + 
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
  
ggsave(paste0("results/3_distance_metrics/heatmaps_all.png"),
       plot = ggall,
       width = 8, height = 20)


# dist_mergeMELT <- 
#   dist_merge %>% 
#   pivot_longer(!c("site", "site1", "site2"), names_to = "variable", values_to = "distance") %>% 
#   dplyr::filter(!grepl('.beta.jne', variable)) %>% 
#   dplyr::filter(!variable %in% c("GstPP.hed", "geodist"))
# 
# dist_mergeMELT$variable <- 
#   factor(dist_mergeMELT$variable,
#          levels = c("seadist", "environment",
#                     "Fst", "D.Jost", 
#                     "Etelinae.beta.jac", "Etelinae.beta.jtu",
#                     "Lutjanidae.beta.jac", "Lutjanidae.beta.jtu" ,    
#                     "Eupercaria/misc.beta.jac", "Eupercaria/misc.beta.jtu",
#                     "Teleostei.beta.jac", "Teleostei.beta.jtu"))
# 
# ggplot(data = dist_mergeMELT, aes(site2, site1, fill = distance))+ 
#   geom_tile()+ 
#   facet_wrap(~variable, scales = "free", ncol = 2) +
#   scale_fill_gradient(low = "#E5F5E0", high = "#31A354") +
#   labs(x = "", y = "") + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
#         panel.grid = element_blank(), 
#         panel.background = element_rect(fill = "white"),
#         plot.background = element_rect(fill = "white"))
# # coord_fixed()
# 
# ggsave(paste0("results/3_distance_metrics/heatmap_ALL.png"),
#        width = 8, height = 20)
# 



# ---- Jaccard partitonning heatmap ----

comm <- "Etelinae"

dist_merge[[paste0(comm, ".ratio")]] <- 
  dist_merge[[paste0(comm, ".beta.jtu")]] /
  dist_merge[[paste0(comm, ".beta.jac")]]
  
gg <- 
  ggplot(data = dist_merge, aes(site2, site1, fill = .data[[paste0(comm, ".ratio")]]))+ 
  geom_tile()+ 
  scale_fill_viridis(direction = -1, option = "mako", begin = 0.05, end = 0.95) +
  labs(fill = "jtu/jac") +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) + 
  coord_fixed()

gg
ggsave(paste0("results/3_distance_metrics/heatmap_SD_", comm, "_ratio_jtu_jne.png"),
       width = 5, height = 4)



