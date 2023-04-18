library(ggrepel)

## ---- parameters ----
level = "site"

# communities delineation
comm_delin = "taxonomic_scale_datasp2"
# comm_delin = "taxonomic_scale_Fishbase"
# comm_delin = "depth_category"
# comm_delin = "phylogenetic_distance"
list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

# parameters
comm = names(list_communities)[1]
metricSD = "richness_site"
metricGD = "Hs"


## ---- load ----
gd_alpha <- read.csv(paste0("results/1_genetic_diversity/gd_table_", level, ".csv"))
sd_alpha <- read.csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))





## ---- SGDCs ----

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
  
  # # filter Seychelles
  # table_alpha <- 
  #   filter(table_alpha, site != "Seychelles")
  
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

gg_grob <- arrangeGrob(grobs = gg_list, ncol=3)
plot(gg_grob)
ggsave(gg_grob, width = 14, height = 5, 
       filename = paste0("results/4_continuity/alpha_SGDCs_all_", level, "_", metricGD, "_", metricSD, "_",  comm_delin, ".png"))


