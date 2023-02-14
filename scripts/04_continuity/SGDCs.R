


## ---- load ----
level = "site"
# level = "station"

list_communities <- readRDS("intermediate/02_species_diversity/List_community.RDS")

gd_global <- read.csv(paste0("results/01_genetic_diversity/gd_table_global_", level, ".csv"))
gd_alpha <- read.csv(paste0("results/01_genetic_diversity/gd_table_", level, ".csv"))
gd_beta <- readRDS(paste0("results/01_genetic_diversity/gd_list_pairwise_", level, ".RDS"))

sd_global <- read.csv("results/02_species_diversity/sd_table_global.csv")
sd_alpha <- read.csv(paste0("results/02_species_diversity/sd_table_", level, ".csv"))
sd_beta <- readRDS(paste0("results/02_species_diversity/sd_list_pairwise_", level, ".RDS"))


# list_communities <- readRDS("intermediate/02_species_diversity/List_community_phylogenetic_scale.RDS")
# sd_beta <- readRDS(paste0("results/02_species_diversity/sd_list_pairwise_", level, "_phylogenetic_scale.RDS"))



## ---- alpha-SGDCs ----

table_alpha <-
  gd_alpha %>% 
  left_join(sd_alpha, by = level)

metricSD = paste0("richness_", level)
metricGD = "Ho"

gg_list <- list()

for (comm in names(list_communities)){
  table_alpha_comm <- 
    table_alpha %>% 
    filter(community == comm)
  
  sgdc_alpha <- summary(lm(table_alpha_comm[[metricGD]] ~ table_alpha_comm[[metricSD]]))
  
  gg_list[[comm]] <- 
    ggplot(table_alpha_comm, aes(.data[[metricSD]], .data[[metricGD]])) +
    geom_point() +
    xlab(paste0("richness_", level, " (", comm, ")")) +
    annotate('text', 
             x=min(table_alpha_comm[[metricSD]]), y=max(table_alpha_comm[[metricGD]]),
             hjust = 0, vjust = 1,
             label=paste0("LM adj.R² = ", sprintf("%.3f", sgdc_alpha$adj.r.squared), "\np = ", sprintf("%.3f", coef(sgdc_alpha)[2,4])))
    
  # ggsave(gg, device = "png", height = 4, width = 6, units = "in",
  #        filename = paste0("results/04_continuity/SGDCs_alpha_all_", level, "_", metricGD, "_", metricSD, "_", comm, ".png"))
  
}

gg_grob <- arrangeGrob(grobs = gg_list, ncol=3)
ggsave(gg_grob, width = 15, height = 8, 
       filename = paste0("results/04_continuity/SGDCs_alpha_all_", level, "_", metricGD, "_", metricSD, ".png"))


## ---- compare Fst, Gst" ----
list_GDbeta <- list()

for (metricGD in names(gd_beta)){
  
  # get distance matrix
  mat_GDbeta <- as.matrix(gd_beta[[metricGD]])

  # order rows alphabetically
  mat_GDbeta <- mat_GDbeta[order(rownames(mat_GDbeta)), order(colnames(mat_GDbeta))]
  
  # pivot longer distance matrix
  melt_GDbeta <- melt.dist(dist = mat_GDbeta, metric = metricGD)
  
  #
  list_GDbeta[[metricGD]] <- melt_GDbeta

}

# merge two distance matrix into one df
GDbeta <- 
  plyr::join_all(list_GDbeta,
            by = c(paste0(level, "1"), paste0(level, "2")))

ggplot(GDbeta, aes(.data[["Fst"]], .data[["GstPP.hed"]])) +
  geom_point()
  






## ---- beta-SGDCs ----

metricSD = "beta.jtu"
metricGD = "Fst"
patt = "all"

gg_list <- list()
comm = names(list_communities)[2]
for (comm in names(list_communities)){
  print(comm)
  
  ## ---- setup beta distance matrix ----
  # get distance matrix
  mat_SDbeta <- as.matrix(sd_beta[[comm]][[metricSD]])
  mat_GDbeta <- as.matrix(gd_beta[[metricGD]])
  
  # keep only sampling sites present in both GD and SD
  mat_SDbeta <- mat_SDbeta[rownames(mat_SDbeta) %in% rownames(mat_GDbeta),
                           colnames(mat_SDbeta) %in% colnames(mat_GDbeta)]
  mat_GDbeta <- mat_GDbeta[rownames(mat_GDbeta) %in% rownames(mat_SDbeta),
                           colnames(mat_GDbeta) %in% colnames(mat_SDbeta)]
  
  # order rows alphabetically
  mat_SDbeta <- mat_SDbeta[order(rownames(mat_SDbeta)), order(colnames(mat_SDbeta))]
  mat_GDbeta <- mat_GDbeta[order(rownames(mat_GDbeta)), order(colnames(mat_GDbeta))]
  
  
  
  ## ---- subset to specific locations ----
  patt = "Seychelles"
  
  mat_SDbeta <- mat_SDbeta[!grepl(patt, rownames(mat_SDbeta)),
                           !grepl(patt, colnames(mat_SDbeta))]
  mat_GDbeta <- mat_GDbeta[!grepl(patt, rownames(mat_GDbeta)),
                           !grepl(patt, colnames(mat_GDbeta))]
  
  
  ## ---- merge and melt beta distance matrix ----
  
  # pivot longer distance matrix
  melt_SDbeta <- melt.dist(dist = mat_SDbeta, metric = metricSD)
  melt_GDbeta <- melt.dist(dist = mat_GDbeta, metric = metricGD)
  
  # merge two distance matrix into one df
  merge_beta <- 
    left_join(melt_SDbeta, melt_GDbeta,
              by = c(paste0(level, "1"), paste0(level, "2")))
  
  # create one column for both locations info
  merge_beta[,level] <- paste0(merge_beta[[paste0(level, "1")]], 
                               "-", 
                               merge_beta[[paste0(level, "2")]])
  
  
  
  ## ---- statistics ----
  # Mantel test
  stat_mantel <- vegan::mantel(as.dist(mat_SDbeta), as.dist(mat_GDbeta))
  print(stat_mantel)
  
  # MRM
  stat_MRM <- MRM(formula = merge_beta[[metricSD]] ~ merge_beta[[metricGD]], nperm = 999)
  print(stat_MRM)
  
  
  ## ---- plot ----
  gg_list[[comm]] <- 
    ggplot(merge_beta, aes(.data[[metricSD]], .data[[metricGD]])) +
    geom_point() +
    xlab(paste0(metricSD, " (", comm, ")")) +
    # geom_text_repel(aes(label = sites), size=3) +
    # ggtitle(paste0(patt)) +
    annotate('text', 
             x=min(merge_beta[[metricSD]]), y=max(merge_beta[[metricGD]]),
             hjust = 0, vjust = 1,
             label=paste0("r Mantel = ", round(stat_mantel$statistic, 4), "\np = ", round(stat_mantel$signif, 4))) 
  
  # ggsave(gg, device = "png", height = 4, width = 6, units = "in",
  #        filename = paste0("results/04_continuity/SGDCs_beta_", patt, "_", level, "_", metricGD, "_", metricSD, "_", comm, ".png"))

}


gg_grob <- arrangeGrob(grobs = gg_list, ncol=3)
ggsave(gg_grob, width = 15, height = 8, 
       filename = paste0("results/04_continuity/SGDCs_beta_noSeychelles_", level, "_", metricGD, "_", metricSD, ".png"))


