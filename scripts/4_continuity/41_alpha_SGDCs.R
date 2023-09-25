library(ggrepel)

# ---- parameters ----
level = "site"

# communtity delineation
# comm_delin = "taxonomy"
# comm_delin = "taxonomy_depth1_crosses45-400m"
# comm_delin = "taxonomy_depth2_contains45-400m"
# comm_delin = "taxonomy_depth3_within45-400m"
comm_delin = "taxonomy_env_reef-associated"

list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

# parameters
comm = names(list_communities)[2]
metricSD = "richness_site"
metricGD = "Hs"


# ---- load ----
gd_alpha <- read_csv(paste0("results/1_genetic_diversity/gd_table_", level, ".csv"))
sd_alpha <- read_csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))





# ---- SGDCs ----

gg_list <- list()

for (comm in names(list_communities)){
  
  # join
  table_alpha <- 
    gd_alpha %>%
    left_join(sd_alpha, by = level) %>% 
    filter(community == comm)
  
  # relevel
  table_alpha[[level]] <- 
    table_alpha[[level]] %>%  
    ordered(levels = unique(data_sites[order(data_sites$order),][[level]]))
  
  # filter Seychelles
  table_alpha <-
    filter(table_alpha, site != "Seychelles")
  
  # LM
  sgdc_alpha <- summary(lm(table_alpha[[metricGD]] ~ table_alpha[[metricSD]]))
  
  # plot
  gg_list[[comm]] <- 
    ggplot(table_alpha, 
           aes(.data[[metricSD]], .data[[metricGD]], 
               color = .data[[level]], label = .data[[level]])) +
    geom_point() +
    geom_text_repel() +
    xlab(paste0("richness_", level, " (", comm, ")")) +
    scale_color_manual(values = color_perso) +
    annotate('text', 
             x=min(table_alpha[[metricSD]]), y=max(table_alpha[[metricGD]]),
             hjust = 0, vjust = 1,
             label=paste0("LM adj.R² = ", sprintf("%.3f", sgdc_alpha$adj.r.squared), "\np = ", sprintf("%.3f", coef(sgdc_alpha)[2,4]))) +
    labs(title = comm) +
    theme(legend.position = "none")
  
  
}

gg_grob <- arrangeGrob(grobs = gg_list, ncol=4)
plot(gg_grob)
ggsave(gg_grob, width = 20, height = 5, 
       filename = paste0("results/4_continuity/alpha_SGDCs_noSeychelles_", level, "_", metricGD, "_", metricSD, "_",  comm_delin, ".png"))


