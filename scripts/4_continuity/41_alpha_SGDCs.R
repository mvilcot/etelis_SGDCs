
# ---- parameters ----
level = "site"

# communtity delineation
comm_delin = "taxonomy"
# comm_delin = "taxonomy_depth1_crosses45-400m"
# comm_delin = "taxonomy_depth2_contains45-400m"
# comm_delin = "taxonomy_depth3_within45-400m"
# comm_delin = "taxonomy_env_reef-associated"

list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

# parameters
comm = names(list_communities)[2]
metricSD = "richness_site"
metricGD = "Hs"


# ---- load ----
# gd_alpha <- read_csv(paste0("results/1_genetic_diversity/gd_table_", level, ".csv"))
# sd_alpha <- read_csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))

alpha_merge <- read_csv(paste0("results/3_distance_metrics/alpha_diversity_", comm_delin, ".csv"))



# ---- SGDCs ----

gg_list <- list()
labs_list <- c("(a)","(b)","(c)","(d)")

for (i in 1:length(list_communities)){
  comm <- names(list_communities)[i]
  lab <- labs_list[i]
  
  # filter Seychelles
  # alpha_sub <-
    # filter(alpha_merge, site != "Seychelles")
  alpha_sub <- alpha_merge
  
  # corr
  # sgdc_alpha <- summary(lm(alpha_sub[[metricGD]] ~ alpha_sub[[metricSD]]))
  sgdc_alpha <- cor.test(alpha_sub[[comm]], alpha_sub[[metricGD]])
  
  # plot
  gg_list[[comm]] <- 
    ggplot(alpha_sub, 
           aes(.data[[comm]], .data[[metricGD]], 
               color = .data[[level]], label = .data[[level]])) +
    geom_point() +
    geom_text_repel() +
    xlab(paste0("species richness (", comm, ")")) +
    scale_color_manual(values = color_perso) +
    annotate('text', 
             x=min(alpha_sub[[comm]]), y=max(alpha_sub[[metricGD]]),
             hjust = 0, vjust = 1,
             # label=paste0("LM adj.R² = ", sprintf("%.3f", sgdc_alpha$adj.r.squared), "\np = ", sprintf("%.3f", coef(sgdc_alpha)[2,4]))) +
             label=paste0("r = ", sprintf("%.3f", sgdc_alpha$estimate), "\np-value = ", sprintf("%.3f", sgdc_alpha$p.value))) +
    labs(tag = lab) + #title = comm
    theme_light() +
    theme(legend.position = "none")
  
  
}

gg_grob <- arrangeGrob(grobs = gg_list, ncol=4)
# plot(gg_grob)
ggsave(gg_grob, width = 20, height = 5, 
       filename = paste0("results/4_continuity/alpha_SGDCs_allsites_", level, "_", metricGD, "_", metricSD, "_",  comm_delin, ".png"))


