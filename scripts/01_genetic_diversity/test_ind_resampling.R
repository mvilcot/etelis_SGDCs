## ---- parameters ----
# filters
filters = "missind_callrate0.70_maf0.05"
level = "site"

# read genlight  
genlight <- 
  read.genlight(filters, level, removeless2ind = TRUE)




## ---- TEST RESAMPLING to 10 ind ----

library(parallel)

N=6

genlight_list <- mclapply(1:10, mc.cores=1, function(i){
  samples_to_keep <- c()
  
  for (pop in levels(genlight$pop)) {
    
    # get sample list from pop
    samples_all <- 
      data_samples %>% 
      dplyr::filter(id %in% genlight$ind.names) %>% 
      dplyr::filter(.[[level]] == pop) %>% 
      pull(id)
    
    # if more than N individuals in a population, sample N individuals names
    if(length(samples_all) > N) { 
      samples_sub <- sample(samples_all, N, replace=F) 
    }
    # if less than N individuals in a population, collect all individual names
    if(length(samples_all) <= N) {
      samples_sub <- samples_all 
    }
    
    # add to sample list to keep
    samples_to_keep <- c(samples_to_keep, samples_sub)
  }
  
  genlight_resampled <- 
    gl.keep.ind(x = genlight, 
                ind.list = samples_to_keep,
                mono.rm = TRUE)
  
})

saveRDS(genlight_list, 
        "intermediate/01_genetic_diversity/TEST/Genlight_resampled_Etelis_coruscans_ordered_missind_callrate0.70_maf0.05.RDS")


for (i in 1:length(genlight_list)){
  cat("0 - LAUNCH", i, "/10")
  
  ## ---- select genlight version ----
  
  i=1
  genlight <- genlight_list[[i]]
  
  
  
  ## ---- convert to other formats ----
  
  # genind 
  genind <- dartR::gl2gi(genlight)
  
  # hierfstat format
  ghierfstat <- hierfstat::genind2hierfstat(genind)
  
  cat("1 - conversion ok", i, "/10")
  
  
  
  ## ---- mean genetic diversity ----
  
  # basic stats
  BS <- hierfstat::basic.stats(genind)
  BS
  BSo <- BS$overall # over all samples
  
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
  
  cat("2 - global metrics ok", i, "/10")
  
  
  
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
  
  
  cat("3 - pairwise metrics ok", i, "/10")
  
  
  ## ---- export ----
  saveRDS(BS, paste0("intermediate/01_genetic_diversity/TEST/basic_stats_", level, "_", i, ".RDS"))
  write.csv(gd_global, paste0("results/01_genetic_diversity/TEST/gd_table_global_", level, "_", i, ".csv"), row.names = F, quote = F)
  write.csv(gd_alpha, paste0("results/01_genetic_diversity/TEST/gd_table_", level, "_", i, ".csv"), row.names = F, quote = F)
  saveRDS(list_gd_beta_pair, paste0("results/01_genetic_diversity/TEST/gd_list_pairwise_", level, "_", i, ".RDS"))
  
  cat("done for", i, "/10")
}



## ---- compare without resampling ----
level = "site"
i=1

gd_global_sub <- read.csv(paste0("results/01_genetic_diversity/TEST/gd_table_global_", level, "_", i, ".csv"))
gd_alpha_sub <- read.csv(paste0("results/01_genetic_diversity/TEST/gd_table_", level, "_", i, ".csv"))
list_gd_beta_pair_sub <- readRDS(paste0("results/01_genetic_diversity/TEST/gd_list_pairwise_", level, "_", i, ".RDS"))

gd_global_all <- read.csv(paste0("results/01_genetic_diversity/gd_table_global_", level, ".csv"))
gd_alpha_all <- read.csv(paste0("results/01_genetic_diversity/gd_table_", level, ".csv"))
list_gd_beta_pair_all <- readRDS(paste0("results/01_genetic_diversity/gd_list_pairwise_", level, ".RDS"))

temp <- t(rbind(gd_global_sub, gd_global_all[,1:13]))
plot(temp[,1], temp[,2])
temp[,1] > temp[,2]

plot(list_gd_beta_pair_sub$Fst, list_gd_beta_pair_all$Fst)
plot(list_gd_beta_pair_sub$GstPP.hed, list_gd_beta_pair_all$GstPP.hed)
plot(list_gd_beta_pair_sub$Fst, list_gd_beta_pair_sub$GstPP.hed)
plot(list_gd_beta_pair_all$Fst, list_gd_beta_pair_all$GstPP.hed)

summary(lm(list_gd_beta_pair_sub$Fst ~ list_gd_beta_pair_all$Fst))
summary(lm(list_gd_beta_pair_sub$GstPP.hed ~ list_gd_beta_pair_all$GstPP.hed))
summary(lm(list_gd_beta_pair_sub$GstPP.hed ~ list_gd_beta_pair_sub$Fst))
summary(lm(list_gd_beta_pair_all$GstPP.hed ~ list_gd_beta_pair_all$Fst))

# Gst is overestimates when low nb and samples, not Fst
list_gd_beta_pair_sub$GstPP.hed > list_gd_beta_pair_all$GstPP.hed
list_gd_beta_pair_sub$Fst > list_gd_beta_pair_all$Fst
