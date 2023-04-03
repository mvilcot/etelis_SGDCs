# read GD and SD genetic diversity
gd_beta <- readRDS("results/TEST/gd_list_pairwise_site_FijiSep.RDS")
sd_beta <- readRDS("results/TEST/sd_list_pairwise_site_depth_category_sepFiji.RDS")

## ---- beta-SGDCs ----
comm_delin = "depth_category"
metricSD = "beta.jtu"
metricGD = "Fst"
patt = "all"

gg_list <- list()
stat_list <- list()
stat_lm <- list()

comm = names(list_communities)[7]
for (comm in names(list_communities)){
  print(comm)
  
  ## -- setup beta distance matrix --
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
  
  
  
  ## -- subset to specific locations --
  patt = "Seychelles"
  
  mat_SDbeta <- mat_SDbeta[!grepl(patt, rownames(mat_SDbeta)),
                           !grepl(patt, colnames(mat_SDbeta))]
  mat_GDbeta <- mat_GDbeta[!grepl(patt, rownames(mat_GDbeta)),
                           !grepl(patt, colnames(mat_GDbeta))]
  
  
  ## -- merge and melt beta distance matrix --
  
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
  
  
  
  ## -- statistics --
  # Mantel test
  stat_mantel <- vegan::mantel(as.dist(mat_SDbeta), as.dist(mat_GDbeta))
  print(stat_mantel)
  
  # LM
  stat_lm[[comm]] <- lm(merge_beta[[metricSD]] ~ merge_beta[[metricGD]])
  
  
  # MRM
  stat_MRM <- MRM(formula = merge_beta[[metricSD]] ~ merge_beta[[metricGD]], nperm = 10000)
  print(stat_MRM)
  stat_list[[comm]] <- stat_MRM$r.squared
  
  ## -- plot --
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
plot(gg_grob)
ggsave(gg_grob, width = 15, height = 8, 
       filename = paste0("results/TEST/SGDCs_beta_noSeychelles_", level, "_", metricGD, "_", metricSD, "_",  comm_delin, ".png"))
