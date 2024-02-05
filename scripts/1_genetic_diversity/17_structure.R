
# ---- read SNPs data set ----
# parameters
filters = "missind1_callrate0.70_maf0.05"
level = "site"
sites = "noSeychelles"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = "Seychelles",
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)


# ---- convert to LEA data format ----

# SNP presence/absence lfmm (package LEA, SilicoDArT)
dartR::gl2geno(genlight,
        outfile = paste0("GenoLEA_Etelis_coruscans_ordered_", sites),
        outpath='intermediate/1_genetic_diversity')


# ---- SNMF ----
kmax = 9

# run snmf
obj.snmf <-
  LEA::snmf(paste0("intermediate/1_genetic_diversity/GenoLEA_Etelis_coruscans_ordered_", sites, ".geno"),
       K = 1:kmax, alpha = 100, project = "new", repetitions = 10, entropy = TRUE)

obj.snmf %>%
  saveRDS(paste0("intermediate/1_genetic_diversity/snmf_Etelis_coruscans_K1-", kmax, "_", sites, ".RDS"))

# plot
pdf(paste0("results/1_genetic_diversity/snmf_Etelis_coruscans_K1-", kmax, "_", sites, "_ggplot2_bis.pdf"),
    height = 4, width = 10)
plot(obj.snmf, cex = 1.2, pch = 19)
for (K in 2:kmax){
  # get the cross-entropy of the 10 runs
  ce = LEA::cross.entropy(obj.snmf, K = K)
  # select the run with the lowest cross-entropy for K = 4
  bestrun = which.min(ce)
  # get best qmatrix
  qmatrix <- LEA::Q(obj.snmf, K = K, run = bestrun)


  ## ---- personalized plot ----
  # create an object with membership probabilities
  probs <- as.data.frame(qmatrix)

  # put probabilities in a tibble with IDS and labels for sites
  probs <-
    as.data.frame(qmatrix) %>%
    cbind(id = genlight@ind.names) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(data_samples[,c("id", "site", "station", "order")], by = "id")

  # melt into long format
  probs_long <-
    probs %>%
    tidyr::pivot_longer(colnames(qmatrix), names_to = "cluster", values_to = "prob")

  # manual relevel of the sampling sites (to avoid alphabetical ordering)
  probs_long$station <- factor(probs_long$station,
                               levels = unique(probs_long[order(probs_long$order),][["station"]]),
                               ordered=TRUE)

  probs_long$site <- factor(probs_long$site,
                            levels = unique(probs_long[order(probs_long$order),][["site"]]),
                            ordered=TRUE)

  # set up custom facet strips
  facetstrips <-
    ggh4x::strip_nested(
      text_x = elem_list_text(size = c(8, 4)),
      by_layer_x = TRUE,
      clip = "off"
    )

  gg <- ggplot(probs_long, aes(factor(id), prob, fill = factor(cluster))) +
    geom_col() +
    ggh4x::facet_nested(~ site,
                 switch = "x",
                 nest_line = element_line(size = 1, lineend = "round"),
                 scales = "free", space = "free", strip = facetstrips) +
    labs(x = "Individuals", y = "membership probability") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_viridis_d() +
    # scale_fill_manual(values = color_perso) +
    theme(
      panel.spacing.x = unit(0.15, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'white', color = 'white'),
      strip.background = element_rect(fill = "white")) +
    labs(fill = "cluster")

  print(gg)

}
dev.off()








# ---- STRUCTURE ----
# install STUCTURE (non GUI version) from here: https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html

## setup genlight ----
genlightLL <- genlight

# add lon/lat (when location by site)
genlightLL@other$latlon <-
  genlightLL@other$ind.metrics %>%
  left_join(shift.lon(data_sites), by = "site") %>%
  dplyr::select(latitude, longitude) %>%
  dplyr::rename(lat = latitude) %>%
  dplyr::rename(lon = longitude)

# subset loci
# nloci <- 500
# genlightLLsub <-
#   gl.subsample.loci(genlightLL, n = nloci, method = 'random', mono.rm = T)


## map lat lon ----
gl.map.interactive(genlightLL)

# change levels to alphabetic, otherwise issue with gl.map.structure
genlightLL$pop <- factor(as.character(genlightLL$pop))


## run structure ----
kmin <- 1
kmax <- 9
nrep <- 5
nloci <- genlightLL$n.loc

structure <- 
  dartR::gl.run.structure(
    genlightLL,
    k.range = kmin:kmax,
    num.k.rep = nrep,
    # pop.prior = "usepopinfo",
    exec = "C:/Users/vilcot/Documents/Structure2.3.4/console/structure.exe",
    plot.out = FALSE)

structure %>% saveRDS(paste0("intermediate/1_genetic_diversity/STRUCTURE_", nloci, "loci_", sites, "_K", kmin, "-", kmax, "_r", nrep, ".RDS"))


## output ---- 
structure <- readRDS(paste0("intermediate/1_genetic_diversity/STRUCTURE_", nloci, "loci_", sites, "_K", kmin, "-", kmax, "_r", nrep, ".RDS"))


# evanno plot
ev <- gl.evanno(structure)

evgg <- # arrange plots
  wrap_plots(ev$plots) +
  patchwork::plot_annotation(tag_levels = 'a', # add tags
                             tag_prefix = '(',
                             tag_suffix = ')',
                             tag_sep = "")

evgg <- # add common x axis
  patchwork::wrap_elements(panel = evgg) +
  labs(tag = "K") + 
  theme(plot.tag = element_text(size = rel(1)),
        plot.tag.position = "bottom")

evgg
ggsave(paste0("results/1_genetic_diversity/STRUCTURE_evanno_", nloci, "loci_", sites, "_K", kmin, "-", kmax, "_r", nrep, ".png"),
       evgg, 
       height = 6, width = 8)


# automatic ancestry plot
k <- 2

# proba matrix
qmat <- gl.plot.structure(structure, K=k, colors_clusters = viridis(k))

# map barplot
gl.map.structure(qmat, genlightLL, K=k, scalex=0.5, scaley=0.2) 



## plot structure ----

### preparing data ----
qmat_plot_large <- 
  as.data.frame(qmat[[1]]) %>% 
  dplyr::select(c("orig.pop", "Label", contains("cluster"))) 
qmat_plot_long <- 
  qmat_plot_large %>% 
  pivot_longer(contains("cluster"), names_to = "cluster", values_to = "ancestry_prop")
qmat_plot_long$orig.pop <- factor(qmat_plot_long$orig.pop, levels = levels(data_sites$site))
qmat_plot_long <- 
  qmat_plot_long %>% 
  dplyr::arrange(orig.pop, ancestry_prop, cluster)



### perso plot ----
gg <- 
  ggplot(qmat_plot_long, aes(Label, ancestry_prop, fill = cluster)) +
  geom_col(width = 1.1) +
  facet_grid(~orig.pop, switch = "x", scales = "free", space = "free") +
  labs(x = "", y = "Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.15)) +
  scale_fill_viridis_d('', labels = c('Cluster 1', 'Cluster 2')) +
  theme(
    panel.spacing.x = unit(0.15, "lines"),
    # axis.text.x = element_text(size=2, angle = 90, vjust = 0.5, hjust=1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.background = element_rect(fill = 'white', color = 'white'),
    strip.background = element_rect(fill = 'white', color = "white"))

gg
# ggsave(paste0("results/1_genetic_diversity/STRUCTURE_barplot_", nloci, "loci_", sites, "_K", k, "_usepopinfo2.png"),
#        gg,
#        height = 4, width = 10)




### perso map ----
# # Proportion of individual from each cluster
# qmat_plot_large <-
#   qmat_plot_large %>% 
#   mutate(CLUSTER = ifelse(cluster1 >= 0.8, "cluster1",
#                           ifelse(cluster2 >= 0.8, "cluster2", "mixed")))
# 
# prop_cluster <-
#   qmat_plot_large %>% 
#   dplyr::rename(site = orig.pop) %>% 
#   group_by(site, CLUSTER) %>% 
#   dplyr::summarise(value = n(), .groups = "keep") %>% 
#   pivot_wider(names_from = CLUSTER, values_from = value)
# prop_cluster$site <- ordered(prop_cluster$site, levels=levels(data_sites$site))
# 
# prop_cluster <-
#   prop_cluster %>% 
#   left_join(data_sites, by = "site")
# 
# prop_cluster <- prop_cluster %>% replace(is.na(.), 0)


# OR Proportion of ancestry in each site
prop_cluster <-
  qmat_plot_large %>% 
  dplyr::rename(site = orig.pop) %>% 
  group_by(site) %>% 
  dplyr::summarise(cluster1 = sum(cluster1), 
                   cluster2 = sum(cluster2), 
                   # cluster3 = sum(cluster3), 
                   # cluster4 = sum(cluster4),
                   # cluster5 = sum(cluster5),
                   .groups = "keep") 
prop_cluster$site <- ordered(prop_cluster$site, levels=levels(data_sites$site))

prop_cluster <-
  prop_cluster %>% 
  left_join(data_sites, by = "site")

# load maps --
# full map
map1 <-
  fortify(maps::map(fill=TRUE, plot=FALSE)) %>% 
  as_tibble()

# add map +360
## >>> SEE IF POSSIBLE TO DO WITH terra::rotate(left = FALSE) #########
map2 <- 
  map1 %>% 
  mutate(long = long + 360,
         group = group + max(map1$group) + 1)

# crop lon & lat extent
map <- 
  rbind(map1, map2) %>% 
  filter(long > 20  & long < 210 & lat <48 & lat > -48)

library(scatterpie)
gg2 <- 
  ggplot() +
  ## Countries
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill="grey20") +
  ## Sites
  geom_scatterpie(data = shift.lon(prop_cluster), 
             aes(x = longitude, y = latitude), #r = number_samples/10
             cols = c("cluster1", "cluster2"), # , "cluster3", "cluster4", "cluster5"
             pie_scale = 1.5, color=NA) + 
  ggrepel::geom_text_repel(data = shift.lon(prop_cluster), show.legend = FALSE,
                           aes(x = longitude, y = latitude, color = site, label = site),
                           hjust=0.5, vjust=0, max.overlaps = 5, nudge_x = -10,
                           bg.color = "grey70", bg.r = 0.02) +
  ## Theme
  scale_color_viridis(discrete = T) +
  scale_fill_viridis('', discrete = T, direction = 1, labels = c('Cluster 1', 'Cluster 2')) +
  theme_minimal() +
  theme(plot.background = element_rect(fill="white", color = "white")) +
  labs(x = "", y = "") +
  coord_equal()
  
gg1
# ggsave(paste0("results/1_genetic_diversity/STRUCTURE_mappiechart_", nloci, "loci_", sites, "_K", k, "_propancestry.png"),
#        gg,
#        height = 8, width = 15)



gg12 <- 
  gg1 / gg2 +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")") 
ggsave(paste0("results/1_genetic_diversity/_S2_STRUCTURE_mappiechart_K", k, "_propancestry2.png"),
       gg12,
       height = 10, width = 10, dpi = 500)





# ---- **** DRAFTS ----


## *** tess3 ----
# # https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html
# library(tess3r)
# 
# genotype <- 
#   read.table(paste0("intermediate/1_genetic_diversity/GenoLEA_Etelis_coruscans_ordered_", sites, ".lfmm"))
# genotype[genotype == 9] <- NA
# 
# coord <- 
#   data_samples %>% 
#   dplyr::filter(id %in% genlight@ind.names) %>% 
#   dplyr::select(Longitude_approx, Latitude_approx) %>% 
#   as.matrix
# 
# # run tess3
# obj.tess3 <- tess3r::tess3(X = genotype, 
#                    coord = coord, K = 1:3,
#                    method = "projected.ls",
#                    ploidy = 2)
# 
# # cross-validation
# plot(obj.tess3, pch = 19, col = "blue",
#      xlab = "Number of ancestral populations",
#      ylab = "Cross-validation score")
# 
# # retrieve tess3 Q matrix 
# q.matrix <- tess3r::qmatrix(obj.tess3, K = 2)
# 
# # STRUCTURE-like barplot for the Q-matrix 
# my.colors <- c("tomato", "orange", "lightblue", "wheat","olivedrab")
# my.palette <- CreatePalette(my.colors, 9)
# barplot(q.matrix, border = NA, space = 0, 
#         main = "Ancestry matrix", 
#         xlab = "Individuals", ylab = "Ancestry proportions", 
#         col.palette = my.palette) -> bp
# axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
# 
# # map plot
# library(rworldmap)
# map.polygon <- getMap(resolution = "low")
# 
# pl <- tess3r::ggtess3Q(q.matrix, coord, map.polygon = map.polygon)
# pl +
#   geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
#   xlim(-180, 180) +
#   ylim(-60, 60) +
#   coord_equal() +
#   geom_point(data = as.data.frame(coord), aes(x = Longitude_approx, y = Latitude_approx), size = 0.2) +
#   xlab("Longitute") +
#   ylab("Latitude") +
#   theme_bw()
# 
# 
## *** TESTS!! ----
# coord = read.table("coordinates.coord")
# pop = rep(1:60, each = 10)
# 
# 
# K = 3 
# Npop = 10
# qpop = matrix(NA, ncol = K, nrow = Npop)
# coord.pop = matrix(NA, ncol = 2, nrow = Npop) 
# for (i in data_sites$site){ 
#   qpop[i,] = apply(qmatrix[pop == i,], 2, mean) 
#   coord.pop[i,] = apply(coord[pop == i,], 2, mean)
#   }
# 
# 
# 
# plot(data_sites, xlab = "longitude", ylab = "latitude", type = "n")
# map(add = T, col = "grey90", fill = TRUE)
# for (i in 1:10){ 
#   add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "", 
#           col = c("orange","violet","lightgreen"))
#   }
# 

## *** dartR test ----
# add lon/lat (when location by site)
genlightLL <- genlight
genlightLL@other$latlon <-
  genlightLL@other$ind.metrics %>%
  left_join(shift.lon(data_sites), by = "site") %>%
  dplyr::select(latitude, longitude) %>%
  dplyr::rename(lat = latitude) %>%
  dplyr::rename(lon = longitude)

nloci <- 1000
genlightLLsub <- 
  gl.subsample.loci(genlightLL, n = nloci, method = 'random', mono.rm = T)

# gl.report.parent.offspring(genlight)
# gl.amova(genlightLLsub)
# genlightLLsub@pop <- factor(rep("ALL", length(genlightLLsub@pop)))
# temp <- gl.LDNe(genlightLLsub, neest.path = "C:/Users/vilcot/Documents/NeEstimator2.1")

beta <- gl.fst.pop(genlight, nboot = 100)

