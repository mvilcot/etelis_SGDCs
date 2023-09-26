
# ---- setup phylogeny ----

# list_species <- 
#   data_species %>% 
#   filter(family %in% "Lutjanidae") %>%
#   pull(species) # keep only species info
# tree_sp  <- fishtree_phylogeny(species = list_species)

tree_lutj  <- fishtree_phylogeny(rank = "Lutjanidae")
tree_lutj$tip.label


# ---- plot parameters ----
# find edges
edge <- which.edge(tree_lutj, "Etelis_coruscans") # target edge
tip <- grep("Etelis_coruscans", tree_lutj$tip.label)

# colour edges
edgecol <- rep('darkgrey', nrow(tree_lutj$edge)) # all default to black
edgecol[edge] <- "chartreuse4"

# thicken edges
edgethick <- rep(1, nrow(tree_lutj$edge))
edgethick[edge] <- 3

# colour tip
tipcol <- rep('darkgrey', length(tree_lutj$tip.label)) # all default to black
tipcol[tip] <- "chartreuse4"

# colour tip
tipfont <- rep(1, length(tree_lutj$tip.label)) # all default to black
tipfont[tip] <- 2



# ---- save plot ----
for (phylotype in c("fan", "tidy", "unrooted")){
  png(paste0("results/0_sampling_design/Phylogeny_lutjanidae_", phylotype, "_raboskytaxo.png"),
      height = 10, width = 10, 
      units = 'in', res = 300)
  
  plot.phylo(tree_lutj,
             type = phylotype,
             align.tip.label = TRUE,
             edge.color = edgecol,
             edge.width = edgethick,
             tip.color = tipcol,
             font = tipfont,
             no.margin = TRUE)
  
  dev.off()
  
}



