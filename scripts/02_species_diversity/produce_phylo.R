## ---- set session ----
# load packages
library(fishtree)
library(readxl)
library(phylotools)


# wrangle phylogenies ----
# sample phylogeny

list_species <- 
  data_species %>% 
  filter(family %in% "Lutjanidae") %>%
  pull(species) # keep only species info

tree_sp  <- fishtree_phylogeny(species = list_species)
# tree_teleo  <- fishtree_phylogeny(rank = "Teleostei")



## ---- plot ----
# find edges
edge <- which.edge(tree_sp, "Etelis_coruscans") # target edge
tip <- grep("Etelis_coruscans", tree_sp$tip.label)

# colour edges
edgecol <- rep('grey', nrow(tree_sp$edge)) # all default to black
edgecol[edge] <- "chartreuse4"

# thicken edges
edgethick <- rep(1, nrow(tree_sp$edge))
edgethick[edge] <- 3

# edge types
edgetype <- rep(1, nrow(tree_sp$edge))
edgetype[edge] <- 1

# colour tip
tipcol <- rep('grey', length(tree_sp$tip.label)) # all default to black
tipcol[tip] <- "chartreuse4"


# jpeg

type = "fan"

png(paste0("results/02_species_diversity/Phylogeny_lutjanidae_", type, ".png"),
    height = 10, width = 10, 
    units = 'in', res = 300)

plot.phylo(tree_sp,
           type = type,
           align.tip.label = TRUE,
           edge.color = edgecol,
           edge.width = edgethick,
           edge.lty = edgetype, 
           tip.color = tipcol,
           no.margin = TRUE)

dev.off()


