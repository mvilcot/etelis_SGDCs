
# ---- read SNPs dataset ----

# parameters
filters = "missind1_callrate0.70_maf0.05"
level = "site"
sites = "noCocos"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL,
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)

genlight



# ---- DAPC pop as prior ----
## run ----
# set.seed(999) # Setting a seed for a consistent result
# dapc <- dapc(genlight)
# saveRDS(dapc, paste0("intermediate/1_genetic_diversity/DAPC_popprior_output_", filters, "_", sites, "_", level, ".RDS"))

# results
dapc <- readRDS(paste0("intermediate/1_genetic_diversity/DAPC_popprior_output_", filters, "_", sites, "_", level, ".RDS"))
# print.dapc(dapc)
# summary.dapc(dapc)
# predict.dapc(dapc)


## alpha score ----
temp <- optim.a.score(dapc)

dapc2 <- dapc(genlight, n.da=50, n.pca=17)



## !!!!!!!!! cross validation ----
#HERE this part doesn't work

mat <- as.data.frame(genlight)
# replace NA by mean
for (i in 1:ncol(mat)){
  mat[,i][is.na(mat[,i])] <- mean(mat[,i], na.rm = TRUE)

}
# mat2 <-
# mat %>% 
#   mutate_if(is.numeric, ~replace_na(., mean(., na.rm = TRUE)))

mat <- tab(genlight, NA.method = "mean")

xval <- adegenet::xvalDapc(x = genlight,
                           grp = pop(genlight), 
                           n.pca.max = 100, 
                           n.rep = 3)
  
  
  # ---- PLOTS ----
  # set up custom facet strips
  facetstrips <- 
    ggh4x::strip_nested(
      text_x = elem_list_text(size = c(8, 4), angle = 90, vjust = 0.5, hjust=1),
      by_layer_x = TRUE,
      clip = "off"
    )
  
  
  ## All sites ----
  dapc <- readRDS(paste0("intermediate/1_genetic_diversity/DAPC_popprior_output_missind1_callrate0.70_maf0.05_noCocos_site_15PCA_2DA.RDS"))
  dapc
  
  ### > scatter plot ----
  # Plot ggplot - adapted from https://github.com/laurabenestan/TD_science_conservation
  
  # get dapc values and metadata
  dapc_geo <- data.frame(id = rownames(dapc$tab), dapc$ind.coord)
  dapc_geo <- 
    data_samples %>% 
    dplyr::right_join(dapc_geo, by = "id") %>% 
    cbind(grp = dapc$grp)
  
  # reorder labels by longitude
  dapc_geo$grp <- factor(dapc_geo$grp, 
                             levels = unique(dapc_geo[order(dapc_geo$order),][[level]]), 
                             ordered=TRUE)
  
  
  # Plot
  gg_scat1 <- ggplot(dapc_geo, aes(x=LD1, y=LD2, color=grp)) +
    geom_point(size = 3, alpha = 0.4) +
    stat_ellipse(level = 0.67) +
    scale_color_manual(values = color_perso) +
    labs(x="DPC1") +
    labs(y="DPC2") +
    # geom_text(label=dapc_geo$id, nudge_x = 0.25, nudge_y = 0.25,
    #           check_overlap = T, show.legend = FALSE)+
    theme_light() +
    theme(legend.position = "none") +
    labs(color = level) +
    labs(tag = "(a)")
  gg_scat1
  
  
  
  ### > compolot ----
  # personalized plot from https://luisdva.github.io/rstats/dapc-plot/
  
  # create an object with membership probabilities
  probs <- as.data.frame(round(dapc$posterior, 4))
  
  # put probabilities in a tibble with IDS and labels for sites
  probs <- 
    tibble::rownames_to_column(probs, var = "id") %>%
    dplyr::left_join(data_samples[,c("id", "site", "station", "order")], by = "id")
  
  # melt into long format
  probs_long <- 
    probs %>% 
    tidyr::pivot_longer(2:(nlevels(dapc$grp)+1), names_to = "cluster", values_to = "prob")
  
  # manual relevel of the sampling sites (to avoid alphabetical ordering)
  probs_long$cluster <- factor(probs_long$cluster, 
                             levels = unique(probs_long[order(probs_long$order),][[level]]), 
                             ordered=TRUE)
  
  probs_long$station <- factor(probs_long$station, 
                               levels = unique(probs_long[order(probs_long$order),][["station"]]), 
                               ordered=TRUE)
  
  probs_long$site <- factor(probs_long$site, 
                               levels = unique(probs_long[order(probs_long$order),][["site"]]), 
                               ordered=TRUE)
  
  
  
  # plot
  gg_compo1 <- ggplot(probs_long, aes(factor(id), prob, fill = factor(cluster))) +
    geom_col(color = "gray", linewidth = 0.01) +
    ggh4x::facet_nested(~ site,
                 switch = "x",
                 nest_line = element_line(linewidth = 1, lineend = "round"),
                 scales = "free", space = "free", strip = facetstrips) +
    labs(x = "Individuals", y = "membership probability") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual(values = color_perso) +
    theme(
      panel.spacing.x = unit(0.15, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'white', color = 'white'),
      strip.background = element_rect(fill = "white")) +
    labs(fill = level) +
    labs(tag = "(c)")
  gg_compo1
  
  
  ## No Seychelles ----
  dapc <- readRDS(paste0("intermediate/1_genetic_diversity/DAPC_popprior_output_missind1_callrate0.70_maf0.05_noCocos_noSeychelles_site.RDS"))
  dapc
  
  ### > scatter plot ----
  # Plot ggplot - adapted from https://github.com/laurabenestan/TD_science_conservation
  
  # get dapc values and metadata
  dapc_geo <- data.frame(id = rownames(dapc$tab), dapc$ind.coord)
  dapc_geo <- 
    data_samples %>% 
    dplyr::right_join(dapc_geo, by = "id") %>% 
    cbind(grp = dapc$grp)
  
  # reorder labels by longitude
  dapc_geo$grp <- factor(dapc_geo$grp, 
                         levels = unique(dapc_geo[order(dapc_geo$order),][[level]]), 
                         ordered=TRUE)
  
  
  # Plot
  gg_scat2 <- ggplot(dapc_geo, aes(x=LD1, y=LD2, color=grp)) +
    geom_point(size = 3, alpha = 0.4) +
    stat_ellipse(level = 0.67) +
    scale_color_manual(values = color_perso) +
    labs(x="DPC1") +
    labs(y="DPC2") +
    # geom_text(label=dapc_geo$id, nudge_x = 0.25, nudge_y = 0.25,
    #           check_overlap = T, show.legend = FALSE)+
    theme_light() +
    theme(legend.position = "none") +
    labs(color = level) +
    labs(tag = "(b)")
  gg_scat2
  
  
  
  ### > compolot ----
  # personalized plot from https://luisdva.github.io/rstats/dapc-plot/
  
  # create an object with membership probabilities
  probs <- as.data.frame(round(dapc$posterior, 4))
  
  # put probabilities in a tibble with IDS and labels for sites
  probs <- 
    tibble::rownames_to_column(probs, var = "id") %>%
    dplyr::left_join(data_samples[,c("id", "site", "station", "order")], by = "id")
  
  # melt into long format
  probs_long <- 
    probs %>% 
    tidyr::pivot_longer(2:(nlevels(dapc$grp)+1), names_to = "cluster", values_to = "prob")
  
  # manual relevel of the sampling sites (to avoid alphabetical ordering)
  probs_long$cluster <- factor(probs_long$cluster, 
                               levels = unique(probs_long[order(probs_long$order),][[level]]), 
                               ordered=TRUE)
  
  probs_long$station <- factor(probs_long$station, 
                               levels = unique(probs_long[order(probs_long$order),][["station"]]), 
                               ordered=TRUE)
  
  probs_long$site <- factor(probs_long$site, 
                            levels = unique(probs_long[order(probs_long$order),][["site"]]), 
                            ordered=TRUE)
  
  
  # plot
  gg_compo2 <- ggplot(probs_long, aes(factor(id), prob, fill = factor(cluster))) +
    geom_col(color = "gray", linewidth = 0.01) +
    ggh4x::facet_nested(~ site,
                        switch = "x",
                        nest_line = element_line(linewidth = 1, lineend = "round"),
                        scales = "free", space = "free", strip = facetstrips) +
    labs(x = "Individuals", y = "membership probability") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual(values = color_perso) +
    theme(
      panel.spacing.x = unit(0.15, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'white', color = 'white'),
      strip.background =element_rect(fill = "white")) +
    theme(legend.position = "none") +
    labs(fill = level) +
    labs(tag = "(d)")
  gg_compo2
  
  
  
  
  
  
  ## Save ----
  gg <- 
    gg_scat1 + gg_scat2 + gg_compo1 + gg_compo2 + 
    plot_layout(heights = c(3, 1), guides = "collect")
  
  ggsave(paste0("results/1_genetic_diversity/DAPC_popprior_perso_", filters, "_", level, ".png"),
         gg, 
         height = 8, width = 14)
  
  
  


# ---- !!!!DAPC with clustering ----
### TO DO!!!!
grp <- find.clusters(genlight, max.n.clust=10)
table(pop(genlight), grp$grp)
dapc1 <- dapc(genlight, grp$grp)
scatter(dapc1, scree.pca = F, legend = T, clabel = F, cleg = 0.8)
compoplot(dapc1, legend = T, show.lab = F, cleg=0.6)

saveRDS(dapc, paste0("intermediate/DAPC_clustering_", filters, "_allsites_5groups.RDS"))








