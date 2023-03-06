

## ---- set parameters ----
# filters
filters <- "missind_callrate0.70_maf0.05"
level <- "site"

# Read genlight  
genlight <- readRDS(paste0("intermediate/01_genetic_diversity/Genlight_Etelis_coruscans_ordered_", filters, ".RDS"))
genlight

# Set population to site/station
genlight@pop <- genlight@other[["ind.metrics"]][[level]]
genlight

# Remove populations with less than two individuals
print(paste("removing", level, ":", names(which(table(genlight@pop) < 2))))
genlight <- gl.drop.pop(genlight, pop.list = names(which(table(genlight@pop) < 2)), recalc = T, mono.rm = T)
genlight



## ---- convert to other formats ----
# Genind 
genind <- dartR::gl2gi(genlight) # same than adegenet::df2genind
genind


## ---- global genetic diversity ----

# basic stats
BS <- hierfstat::basic.stats(genind) # same result
BS
BSo <- BS$overall # over all samples

# # Jost Additive framework
# Hs = BSo[["Hs"]]
# Ht = BSo[["Ht"]]
# Hst = (Ht-Hs)/(1-Hs)
# Dst <- Ht - Hs
# 
# # Jost Multiplicative framework Jst=Jt*Js
# Js = 1/(1-Hs)
# Jt = 1/(1-Ht)
# Jst = Jt/Js # or 1/(1-Hst)

# Jost D
JostD <- D_Jost(genind)

# # Hedrick Gst"
# GstPP.hed <- Gst_Hedrick(genind)




## ---- beta gd pairwise ----

# # Fst
# Fst_pair <- genet.dist(genind, method = "WC84")
# 
# # G"st
# GstPP.hed_pair <- pairwise_Gst_Hedrick(genind)

# Jost D
JostD_pair <- pairwise_D(genind)




## ---- export ----
# put into table
gd_global <- read.csv(paste0("results/01_genetic_diversity/gd_table_global_", level, ".csv"))

gd_global <- 
  gd_global %>% 
  cbind(JostD = JostD$global.het) %>%
  as_tibble()

# put into list
list_gd_beta_pair <- readRDS(paste0("results/01_genetic_diversity/gd_list_pairwise_", level, ".RDS"))

list_gd_beta_pair[["JostD"]] <- 
  JostD_pair

# save
saveRDS(BS, "intermediate/01_genetic_diversity/basic_stats.RDS")
write.csv(gd_global, paste0("results/01_genetic_diversity/gd_table_global_", level, "_withJostD.csv"), row.names = F, quote = F)
saveRDS(list_gd_beta_pair, paste0("results/01_genetic_diversity/gd_list_pairwise_", level, "_withJostD.RDS"))



