
## ---- parameters ----
# filters
filters = "missind_callrate0.70_maf0.05"
level = "site"

# read genlight  
genlight <- 
  read.genlight(filters, level, removeless2ind = TRUE)

########### SPECIFIC TO FIJI ###############
data_samples <- read.csv("scripts/TESTS/metadata_samples_Fiji.csv")
data_sites <- read.csv("scripts/TESTS/metadata_sites_Fiji.csv")

data_samples <- 
  data_samples %>% 
  left_join(data_sites, by = c("station", "site"))


# replace metadata
genlight@other[["ind.metrics"]] <- 
  data_samplesFiji %>% 
  dplyr::filter(id %in% genlight@other$ind.metrics$id)

# replace pop info
genlight@pop <- 
  genlight@other[["ind.metrics"]][[level]] %>% 
  as.factor()


## ---- convert to other formats ----

# genind 
genind <- dartR::gl2gi(genlight)

# hierfstat format
ghierfstat <- hierfstat::genind2hierfstat(genind)



## ---- mean genetic diversity ----

# basic stats
BS <- hierfstat::basic.stats(genind)
BSo <- BS$overall

# Jost additive framework Dst = Ht - Hs
Hs = BSo[["Hs"]]
Ht = BSo[["Ht"]]
Hst = (Ht-Hs)/(1-Hs)
Dst <- Ht - Hs

# Jost multiplicative framework Jst = Jt / Js
Js = 1/(1-Hs)
Jt = 1/(1-Ht)
Jst = Jt/Js # or 1/(1-Hst)

# Hedrick Gst" (Meirmans and Hedrick 2011)
GstPP.hed <- Gst_Hedrick(genind)

# Jost D (Jost 2008)
JostD <- D_Jost(genind)

# Fst population-specific (Weir and Goudet 2017)
popFst <- betas(ghierfstat, nboot=1000)

# create global table
gd_global <-
  cbind(Ho = BSo["Ho"], Hs, Ht, Hst, Dst,
        Js, Jt, Jst,
        Fst = BSo["Fst"], 
        Fis = BSo["Fis"],
        Dstp = BSo["Dstp"], 
        Dest = BSo["Dest"],
        Gstpp.hed = GstPP.hed$global,
        D.Jost = JostD$global.het,
        popFst.WG = popFst$betaW) %>%
  as_tibble()




## ---- alpha gd by location ----
gd_alpha <-
  data.frame(Hs = colMeans(BS$Hs, na.rm = T),
             Ho = colMeans(BS$Ho, na.rm = T)) %>% 
  rownames_to_column(level) %>% 
  left_join(
    data.frame(popFst.WG = popFst$betaiovl) %>% 
      rownames_to_column(level), 
    by = level)




## ---- beta gd pairwise ----

# Fst
Fst_pair <- genet.dist(genind, method = "WC84")

# Hedrick G"st
GstPP.hed_pair <- pairwise_Gst_Hedrick(genind)

# Jost D
JostD_pair <- pairwise_D(genind)

# put into list
list_gd_beta_pair <-
  list(Fst = Fst_pair, 
       GstPP.hed = GstPP.hed_pair,
       D.Jost = JostD_pair)


## ---- export ----
saveRDS(BS, paste0("intermediate/TEST/basic_stats_", level, "_FijiSep.RDS"))
write.csv(gd_global, paste0("results/TEST/gd_table_global_", level, "_FijiSep.csv"), row.names = F, quote = F)
write.csv(gd_alpha, paste0("results/TEST/gd_table_", level, "_FijiSep.csv"), row.names = F, quote = F)
saveRDS(list_gd_beta_pair, paste0("results/TEST/gd_list_pairwise_", level, "_FijiSep.RDS"))


