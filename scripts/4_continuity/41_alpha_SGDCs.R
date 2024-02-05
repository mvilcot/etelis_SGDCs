
# ---- parameters ----
level = "site"

# communtity delineation
comm_delin_list <-
  c("taxonomy",
    "taxonomy_depth1_crosses45-400m",
    "taxonomy_depth2_contains45-400m",
    "taxonomy_depth3_within45-400m",
    "taxonomy_env_reef-associated")
comm_delin <- comm_delin_list[1]

list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))
names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")

# parameters
metricSD = "richness_site"
metricGD = "Hs"




# ---- table ----
i=1
metricGD = "Hs"
SGDC_alpha <- tibble(community_delineation = NA,
                     locations = NA,
                     taxonomic_scale = NA,
                     metricGD = NA,
                     metricSD = NA,
                     r = NA,
                     pval = NA
)

for(comm_delin in comm_delin_list) {
  
  alpha_merge <- read_csv(paste0("results/3_distance_metrics/alpha_diversity_", comm_delin, ".csv"))
  
  for (locations in c("all sites", "without Seychelles")){
    
    if(locations == "without Seychelles"){
      alpha_sub <- filter(alpha_merge, site != "Seychelles")
    }
    alpha_sub <- alpha_merge

    
    for (metricSD in names_communities){
    
        # corr
        stat_pearson <- cor.test(alpha_sub[[metricGD]], alpha_sub[[metricSD]])
        
        SGDC_alpha[i,]$community_delineation <- comm_delin
        SGDC_alpha[i,]$locations <- locations
        SGDC_alpha[i,]$taxonomic_scale <- metricSD
        SGDC_alpha[i,]$metricGD <- metricGD
        SGDC_alpha[i,]$metricSD <- "species richness"
        SGDC_alpha[i,]$r <- stat_pearson$estimate
        SGDC_alpha[i,]$pval <- stat_pearson$p.value
        
        cat(i, "\n")
        i=i+1
    
    }
  }
}


SGDC_alpha$signif <- ifelse(SGDC_alpha$pval < 0.001, "***", 
                            ifelse(SGDC_alpha$pval < 0.01, "**", 
                                   ifelse(SGDC_alpha$pval < 0.05, "*", "NS")))

SGDC_alpha %>%
  write_csv("results/4_continuity/alpha_SGDCs_table.csv")




# ---- plot ----
comm_delin <- comm_delin_list[1]
alpha_merge <- read_csv(paste0("results/3_distance_metrics/alpha_diversity_", comm_delin, ".csv"))
alpha_merge$site <- factor(alpha_merge$site,
                           level = unique(alpha_merge$site))

gg_list <- list()
for (comm in names_communities){
  # filter Seychelles
  alpha_sub <- filter(alpha_merge, site != "Seychelles")
  LABELS_sub <- LABELS[-1]
  # alpha_sub <- alpha_merge
  # LABELS_sub <- LABELS
  
  # corr
  # stat_pearson <- summary(lm(alpha_sub[[metricGD]] ~ alpha_sub[[metricSD]]))
  stat_pearson <- cor.test(alpha_sub[[comm]], alpha_sub[[metricGD]])
  
  # plot
  gg_list[[comm]] <-
    ggplot(alpha_sub, 
           aes(.data[[comm]], .data[[metricGD]], 
               color = .data[[level]], label = .data[[level]])) +
    geom_point(size = 3, alpha = 0.7) +
    # geom_text_repel(bg.color = "grey70", bg.r = 0.02, max.overlaps = 0.5, show.legend = F) +
    xlab(paste0("species richness (", comm, ")")) +
    scale_color_manual('', values = color_perso, labels = LABELS_sub) +
    annotate('text', 
             x=min(alpha_sub[[comm]]), y=1.01*max(alpha_sub[[metricGD]]),
             hjust = 0, vjust = 1,
             label=paste0("r Pearson = ", round(stat_pearson$estimate, 3), "\np = ", round(stat_pearson$p.value, 3))) +
    theme_light()
    # theme(legend.position = "none")
  
  
}

### {FIGURE S6} ####
gg_grob <-
  patchwork::wrap_plots(gg_list, tag_level = "new", guides = "collect") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
plot(gg_grob)
ggsave(gg_grob, width = 9, height = 7, dpi = 500,
       filename = paste0("results/4_continuity/_S6_alpha_SGDCs_noSeychelles_plot_", comm_delin, ".png"))


