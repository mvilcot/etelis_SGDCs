# library(fields)
# library(RColorBrewer)
# library(mapplots)



# ---- read SNPs dataset ----

# parameters
filters = "missind1_callrate0.70_maf0.05"
level = "site"
sites = "all"
nrun = 9

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL,
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)


# ---- convert to LEA data format ----

# SNP presence/absence lfmm (package LEA, SilicoDArT)
dartR::gl2geno(genlight, 
        outfile = paste0("GenoLEA_Etelis_coruscans_ordered_", sites),
        outpath='intermediate/1_genetic_diversity')


# ---- snmf ----
library(LEA)

# run snmf
obj.snmf <-
  LEA::snmf(paste0("intermediate/1_genetic_diversity/GenoLEA_Etelis_coruscans_ordered_", sites, ".geno"), 
       K = 1:nrun, alpha = 100, project = "new", repetitions = 10, entropy = TRUE)

obj.snmf %>% 
  saveRDS(paste0("intermediate/1_genetic_diversity/snmf_Etelis_coruscans_K1-", nrun, "_", sites, ".RDS"))

# plot 
pdf(paste0("results/1_genetic_diversity/snmf_Etelis_coruscans_K1-", nrun, "_", sites, "_ggplot2.pdf"),
    height = 4, width = 10)
plot(obj.snmf, cex = 1.2, pch = 19)
for (K in 2:nrun){
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






## TESTS!!
# coord = read.table("coordinates.coord") 
# pop = rep(1:60, each = 10)
# 
# 
# K = 3 
# Npop = 10
# qpop = matrix(NA, ncol = K, nrow = Npop)
# coord.pop = matrix(NA, ncol = 2, nrow = Npop) 
# for (i in rownames(coord_site)){ 
#   qpop[i,] = apply(qmatrix[pop == i,], 2, mean) 
#   coord.pop[i,] = apply(coord[pop == i,], 2, mean)
#   }
# 
# 
# 
# plot(coord_site, xlab = "longitude", ylab = "latitude", type = "n")
# map(add = T, col = "grey90", fill = TRUE)
# for (i in 1:10){ 
#   add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "", 
#           col = c("orange","violet","lightgreen"))
#   }
# 



# ---- tess3 ----
library(tess3r)

genotype <- 
  read.table(paste0("intermediate/1_genetic_diversity/GenoLEA_Etelis_coruscans_ordered_", sites, ".lfmm"))
genotype[genotype == 9] <- NA

coord <- 
  data_samples %>% 
  dplyr::filter(id %in% genlight@ind.names) %>% 
  dplyr::select(Longitude_approx, Latitude_approx) %>% 
  as.matrix

# run tess3
obj.tess3 <- tess3r::tess3(X = genotype, 
                   coord = coord, K = 1:3,
                   method = "projected.ls",
                   ploidy = 2)

# cross-validation
plot(obj.tess3, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# retrieve tess3 Q matrix 
q.matrix <- tess3r::qmatrix(obj.tess3, K = 2)

# STRUCTURE-like barplot for the Q-matrix 
my.colors <- c("tomato", "orange", "lightblue", "wheat","olivedrab")
my.palette <- CreatePalette(my.colors, 9)
barplot(q.matrix, border = NA, space = 0, 
        main = "Ancestry matrix", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 

# map plot
library(rworldmap)
map.polygon <- getMap(resolution = "low")

pl <- tess3r::ggtess3Q(q.matrix, coord, map.polygon = map.polygon)
pl +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  xlim(-180, 180) +
  ylim(-60, 60) +
  coord_equal() +
  geom_point(data = as.data.frame(coord), aes(x = Longitude_approx, y = Latitude_approx), size = 0.2) +
  xlab("Longitute") +
  ylab("Latitude") +
  theme_bw()

# ---- STRUCTURE ----
library(strataG)
structure <- dartR::gl.run.structure(
  genlight,
  k.range = 2:5,
  num.k.rep = 3,
  # exec = "C:/Users/User/Downloads/structure_windows_console/console/structure.exe",
  plot.out = TRUE,
  save2tmp = FALSE,
  verbose = NULL
)
