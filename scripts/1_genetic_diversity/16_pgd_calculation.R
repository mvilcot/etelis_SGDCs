# caluclate phylogenetic metrics on intra-species phylogenies


## ---- load data ----
# read genlight
filters = "missind_callrate0.70_maf0.05"
level = "site"

genlight <- 
  read.genlight(filters, level, removeless2ind = FALSE)

# convert to genind
genind <- dartR::gl2gi(genlight) # same than adegenet::df2genind
genind


## ---- calculate tree ----
# This script attempts to create a "phylogenetic" tree out of genotype data 
# using pairwise distances of loci between individuals and naive clustering 
# methods

# convert to loci format
data_loci <- 
  pegas::genind2loci(genind)

# generate individual-individual distance matrices
dist_loci <-
  ape::dist.gene(data_loci,
            method = "percentage",
            pairwise.deletion = TRUE, # remove missing loci
            variance = FALSE)

# calculate "phylogenetic" tree using ape neighbour-joining method
tree_gd <- ape::nj(dist_loci)

# export
tree_gd %>% 
  saveRDS("intermediate/1_genetic_diversity/gd_phylo_global.RDS")



## ---- plot tree ----

tree_gd <- readRDS("intermediate/1_genetic_diversity/gd_phylo_global.RDS")

# add metadata
tree_meta <-
  data.frame(id = tree_gd$tip.label) %>%  # get samples in tree
  dplyr::left_join(data_samples, by = "id")      # join metadata

# order sites
tree_meta[[level]] <- factor(tree_meta[[level]],
                                 levels = unique(tree_meta[order(tree_meta$order),][[level]]),
                                 ordered=TRUE)

# plot tree
ggtree(tree_gd) %<+%
  tree_meta +
  geom_tiplab(size=1, offset = 0.001) +
  geom_tippoint(aes(color = site), size=3, alpha=.75) +
  scale_color_viridis_d(paste0(level, "s")) +
  theme_tree2()

# save
ggsave("results/1_genetic_diversity/Tree_nj_genetic_distance.png",
       height = 20, width = 8)


# sites <- factor(tree_meta[[level]], 
#                 levels = unique(tree_meta[order(tree_meta$order),][[level]]), 
#                 ordered=TRUE)
# 
# color_site <-
#   data.frame(site = levels(sites),
#              colour = viridis(length(levels(sites))))
# 
# color_samples <-
#   tree_meta %>% 
#   left_join(color_site, by = level)
# 
# color_samples[[level]] <- factor(color_samples[[level]], 
#                                  levels = unique(color_samples[order(tree_meta$order),][[level]]), 
#                                  ordered=TRUE)
# 
# plot.phylo(tree_gd,
#            tip.color = color_samples$colour)
# legend("topright", title=level,
#        legend = color_site$site,
#        fill = color_site$colour, cex=0.2)




## ---- calculate pgd ----

# create psuedo-community matrix (all ind present)
com_mat <-
  matrix(rep(1,length(tree_gd$tip.label)),
         nrow = 1)

colnames(com_mat) <- tree_gd$tip.label
rownames(com_mat) <- c("global")

# midpoint root the tree
tree_gd <- phytools::midpoint.root(tree_gd)

# pd, mpd, vpd
phylo_metrics <-
  data.frame(gd_pd   = as.numeric(picante::pd(com_mat,tree_gd)[1]),
             gd_mpd  = mean(dist_loci),
             gd_vpd  = var(dist_loci))

## ---- export ----
phylo_metrics %>% 
  write_csv("intermediate/1_genetic_diversity/gd_phylo_metrics.csv")

