
##########!!!!! NEED TO AUMOMATIZE PARAMETERS WITH FUNCTION !!!!! ##############

## ---- read SNPs dataset ----

# set parameters
filters <- "missind_callrate0.70_maf0.05"

# sites <- "allsites"
sites <- "noSeychelles"
# sites <- "Hawaii"
# sites <- "NCaledonia"
# sites <- "WAustralia"

# level <- "site"
level <- "station"

# read genlight
genlight <- readRDS(paste0("intermediate/01_genetic_diversity/Genlight_Etelis_coruscans_ordered_", filters, ".RDS"))
genlight

# drop or keep pop
# genlight <- gl.drop.pop(genlight, pop.list = c("Seychelles"), recalc = T, mono.rm = T)
# genlight <- gl.keep.pop(genlight, pop.list = c("Hawaii"), recalc = T, mono.rm = T)
# genlight <- gl.keep.pop(genlight, pop.list = c("New_Caledonia"), recalc = T, mono.rm = T)
# genlight <- gl.keep.pop(genlight, pop.list = c("W_Australia"), recalc = T, mono.rm = T)
genlight@pop <-
  genlight@other[["ind.metrics"]][[level]] %>%
  droplevels()


## ---- DAPC ----

# DAPC - populations as prior
# set.seed(999) # Setting a seed for a consistent result
# dapc <- dapc(genlight)
# saveRDS(dapc, paste0("intermediate/01_genetic_diversity/DAPC_popprior_output_", filters, "_", sites, "_", level, ".RDS"))

# results
dapc <- readRDS(paste0("intermediate/01_genetic_diversity/DAPC_popprior_output_", filters, "_", sites, "_", level, ".RDS"))
# print.dapc(dapc)
# summary.dapc(dapc)
# predict.dapc(dapc)


## ---- scatter plot ----
# png(paste0("results/01_genetic_diversity/DAPC_popprior_scatter_", filters, "_", sites, "_", level, ".png"),
#     height = 12, width = 12,
#     units = 'in', res = 300)
# scatter(dapc, scree.pca = F, legend = T, clabel = F, cleg = 0.8)
# dev.off()


# Plot ggplot - adapted from https://github.com/laurabenestan/TD_science_conservation

# get dapc values and metadata
# dapc <- readRDS(paste0("intermediate/01_genetic_diversity/DAPC_popprior_output_", filters, "_", sites, ".RDS"))
dapc_geo <- data.frame(id = rownames(dapc$tab), dapc$ind.coord)
dapc_geo <- 
  data_samples %>% 
  right_join(dapc_geo, by = "id") %>% 
  cbind(grp = dapc$grp)

# reorder labels by longitude
dapc_geo$grp <- factor(dapc_geo$grp, 
                           levels = unique(dapc_geo[order(dapc_geo$order),][[level]]), 
                           ordered=TRUE)


# Plot
gg1 <- ggplot(dapc_geo, aes(x=LD1, y=LD2, color=grp)) +
  geom_point(size = 3, alpha = 0.4) +
  stat_ellipse(level = 0.67) +
  # scale_color_manual(values=wes_palette('Zissou1', nlevels(dapc_geo$station), 'continuous')) +
  # scale_color_paletteer_d("ggthemes::Hue_Circle") +
  scale_color_viridis_d() +
  labs(x="DPC1") +
  labs(y="DPC2") +
  # geom_text(label=dapc_geo$id,
  #   nudge_x = 0.25, nudge_y = 0.25,
  #   check_overlap = T, show.legend = FALSE)+
  # coord_fixed(ratio = 1) +
  theme_classic() +
  labs(tag = "A")
# gg1



## ---- compolot plot ----

# png(paste0("results/01_genetic_diversity/DAPC_popprior_compoplot_", filters, "_", sites, "_", level, ".png"),
#     height = 8, width = 20,
#     units = 'in', res = 100)
# compoplot(dapc, legend = T, show.lab = F, cleg=0.6)
# dev.off()

# personalized plot from https://luisdva.github.io/rstats/dapc-plot/

# create an object with membership probabilities
probs <- as.data.frame(round(dapc$posterior, 4))

# put probabilities in a tibble with IDS and labels for sites
probs <- 
  tibble::rownames_to_column(probs, var = "id") %>%
  left_join(data_samples[,c("id", "site", "station", "order")], by = "id")

# melt into long format
probs_long <- 
  probs %>% 
  pivot_longer(2:(nlevels(dapc$grp)+1), names_to = "cluster", values_to = "prob")

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


# set up custom facet strips
facetstrips <- 
  strip_nested(
    text_x = elem_list_text(size = c(8, 4)),
    by_layer_x = TRUE,
    clip = "off"
  )

# plot
gg2 <- ggplot(probs_long, aes(factor(id), prob, fill = factor(cluster))) +
  geom_col(color = "gray", size = 0.01) +
  facet_nested(~ site + station,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips) +
  labs(x = "Individuals", y = "membership probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  # scale_fill_paletteer_d("ggthemes::Hue_Circle") +
  # scale_fill_hue() +
  scale_fill_viridis_d() +
  theme(
    panel.spacing.x = unit(0.15, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = 'white', color = 'white'),
    strip.background =element_rect(fill = "white")) +
  labs(tag = "B")
# gg2



# save
gg1 / gg2 + plot_layout(heights = c(3, 1))
ggsave(paste0("results/01_genetic_diversity/DAPC_popprior_perso_", filters, "_", sites, "_", level, ".png"),
              height = 12, width = 12)




# ## ----------- TO DO!! - DAPC with clustering ---------------------------------
# grp <- find.clusters(genlight, max.n.clust=20)
# table(pop(genlight), grp$grp)
# dapc1 <- dapc(genlight, grp$grp)
# scatter(dapc1, scree.pca = F, legend = T, clabel = F, cleg = 0.8)
# compoplot(dapc1, legend = T, show.lab = F, cleg=0.6)
# 
# saveRDS(dapc, paste0("intermediate/DAPC_clustering_", filters, "_allsites_5groups.RDS"))
# 
# 
# 
# 
# ## ----------- TO DO!! - PCA --------------------------------------------------
# 
# PCA <- glPca(genlight, center = TRUE, scale = FALSE, nf = 10, loadings = TRUE, 
#              alleleAsUnit = FALSE, useC = TRUE, parallel = require("parallel"),
#              n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)
# 
# DataPlot <- as.data.frame(PCA$scores)
# DataPlot$pop <- genlight$pop
# ggplot(DataPlot) +
#   geom_point(aes(PC1, PC2, color=pop)) +
#   theme_classic()
# 
# 



