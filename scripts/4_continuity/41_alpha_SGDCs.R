
# ---- parameters ----
level = "site"

# communtity delineation
comm_delin <- "taxonomy"

list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))
names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")

# parameters
metricSD = "richness_site"
metricGD = "Hs"



# ---- alpha-SGDC ----
## table ----
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


SGDC_alpha$signif <- ifelse(SGDC_alpha$pval < 0.001, "***", 
                            ifelse(SGDC_alpha$pval < 0.01, "**", 
                                   ifelse(SGDC_alpha$pval < 0.05, "*", "NS")))

SGDC_alpha %>%
  write_csv("results/4_continuity/alpha_SGDCs_table.csv")




## plot ----
alpha_merge <- read_csv(paste0("results/3_distance_metrics/alpha_diversity_", comm_delin, ".csv"))
alpha_merge$site <- factor(alpha_merge$site,
                           level = unique(alpha_merge$site))

gg_list <- list()
for (comm in names_communities){
  # filter Seychelles
  alpha_sub <- filter(alpha_merge, site != "Seychelles")
  LABELS_sub <- LABELS[-1]
  
  # corr
  stat_pearson <- cor.test(alpha_sub[[comm]], alpha_sub[[metricGD]])
  
  # plot
  gg_list[[comm]] <-
    ggplot(alpha_sub, 
           aes(.data[[comm]], .data[[metricGD]], 
               color = .data[[level]], label = .data[[level]])) +
    geom_point(size = 3, alpha = 0.7) +
    xlab(paste0("species richness (", comm, ")")) +
    scale_color_manual('', values = color_perso, labels = LABELS_sub) +
    annotate('text', 
             x=min(alpha_sub[[comm]]), y=1.01*max(alpha_sub[[metricGD]]),
             hjust = 0, vjust = 1,
             label=paste0("r = ", round(stat_pearson$estimate, 3), "\np = ", round(stat_pearson$p.value, 3))) +
    theme_light()
  
}

### {FIGURE S7} ####
gg_grob <-
  patchwork::wrap_plots(gg_list, tag_level = "new", guides = "collect") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
plot(gg_grob)
ggsave(gg_grob, width = 9, height = 7, dpi = 500,
       filename = paste0("results/4_continuity/_S7_alpha_SGDCs_noSeychelles_plot_", comm_delin, ".png"))









# Richness ~ longitude ----

## distance to coral triangle ----

# coordinates of the center of the Coral Triangle
center_CT <- 
  data.frame(longitude=133.679826, 
             latitude=-1.307436)

# compute distance from the different stations to the CT (in meters)
dist_to_CT <-
  data_sites %>% 
  cbind(data.frame(distCT = pointDistance(data_sites[,2:3], center_CT, lonlat=TRUE)))



## setup data ----
level = "site"
alpha_merge <- read_csv(paste0("results/3_distance_metrics/alpha_diversity_", comm_delin, ".csv"))

df <-
  dist_to_CT %>% 
  left_join(alpha_merge, by = "site") %>% 
  dplyr::select(-c(Ho, popFst.WG, number_samples))

df$site <- factor(df$site, levels = unique(df$site))
df <- shift.lon(df)

df_long <- 
  df %>% 
  pivot_longer(cols = -c(site, longitude, latitude, distCT), names_to = "variable") %>% 
  mutate(variable = gsub("Hs", "Etelis coruscans (Hs)", variable)) %>% 
  mutate(variable = factor(variable, 
                           levels = c(names_communities, "Etelis coruscans (Hs)"))) %>% 
  cbind(locations = 'All sites')

df_long <-
  df_long %>%
  filter(site != "Seychelles") %>%
  mutate(locations = "Without Seychelles") %>%
  rbind(df_long)




## plot ----
gglist <- list()
for(var in levels(df_long$variable)){
  df_sub <- 
    df_long %>% 
    dplyr::filter(variable == var)
  
  df_subALL <- 
    df_sub %>% 
    dplyr::filter(locations == 'All sites')
  
  df_subNS <- 
    df_sub %>% 
    dplyr::filter(locations == "Without Seychelles")
  
  LMall <- cor.test(df_subALL$value, df_subALL$longitude) # distCT
  LMns <- cor.test(df_subNS$value, df_subNS$longitude) # distCT
  
  gglist[[var]] <-
    df_sub %>% 
    ggplot(aes(x=longitude, y=value, color=site, linetype = locations)) + 
    facet_wrap( ~ variable, scales ="free_y", ncol = 2) +
    geom_smooth(method = "lm", color = "grey40", se = F) +
    geom_point() +
    scale_color_manual('', values = color_perso, labels = LABELS) +
    annotate('text', cex = 2.5,
             x = min(df_sub$longitude), y = 0.9*min(df_sub$value),
             hjust = 0, vjust = 0,
             label = paste0("All sites: r = ", round(LMall$estimate, 4), "; p = ", round(LMall$p.value,3), "\n",
                            "Without Seychelles: r = ", round(LMns$estimate, 4), "; p = ", round(LMns$p.value, 3))) +
    ylab('α-diversity') +
    theme_light() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank())
}

gg_grob <-
  patchwork::wrap_plots(gglist, guides = 'collect', ncol = 2) +
  xlab('Distance to Coral Triangle')

plot(gg_grob)


### {FIGURE S1} ####
ggsave(paste0("results/4_continuity/_S1_alpha_sd_richness_longitude_site.png"),
       width = 9, height = 9, dpi = 500)
ggsave(paste0("results/4_continuity/_S1_alpha_sd_richness_longitude_site.pdf"),
       width = 9, height = 9)

