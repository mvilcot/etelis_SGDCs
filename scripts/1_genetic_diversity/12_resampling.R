
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


# ---- Method 1 : pop by pop, Hs ~ sampling size ----
## ---- downsample genlight ----
genlight_list <- list()

for (pop in levels(genlight$pop)) {
  genlight_temp <- list()
  
  # get samples name for the pop considered
  samples_all <- 
    data_samples %>% 
    dplyr::filter(id %in% genlight$ind.names) %>% 
    dplyr::filter(.data[[level]] == pop) %>% 
    dplyr::pull(id)
  
  for (N in 1:length(samples_all)){
    
    # if more than N individuals in a pop, randomly sample N individuals names
    if(length(samples_all) > N) { 
      samples_to_keep <- sample(samples_all, N, replace=F) 
    }
    
    # if equal to N individuals in a pop, collect all individual names
    if(length(samples_all) <= N) {
      samples_to_keep <- samples_all 
    }
    
    # resample genlight
    genlight_temp[[N]] <- 
      dartR::gl.keep.ind(x = genlight, 
                         ind.list = samples_to_keep,
                         mono.rm = TRUE)
  }

  genlight_list[[pop]] <- genlight_temp

}

genlight_list %>% 
  saveRDS("intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered_missind_callrate0.70_maf0.05_different_sampling_size_POPBYPOP.RDS")



## ---- iteratively compute Hs ----

het_list <- list()
for (pop in levels(genlight$pop)) {
  het_temp <- list()
  genlight_temp <- genlight_list[[pop]]
  
  for (N in 1:length(genlight_temp)){
    het_temp[[N]] <- 
      dartR::gl.report.heterozygosity(genlight_temp[[N]], plot.out = F) %>% 
      mutate(N = N)
  }
  het_list[[pop]] <- 
    het_temp %>% 
    reduce(rbind)
}

het_list %>% saveRDS("intermediate/1_genetic_diversity/Alpha_diversity_dartR_different_sampling_size_POPBYPOP.RDS")


het_df <-
  het_list %>% 
  reduce(rbind)

het_df$site <- factor(x = het_df$pop, 
                      levels = levels(genlight@pop), 
                      ordered = TRUE)

## ---- plot saturation curve ----
ggplot(het_df, aes(N, He, color=site)) +
  scale_color_manual(values = color_perso) +
  geom_line(size = 0.8, alpha = 0.5) +
  xlim(0, 30)
ggsave(paste0("results/1_genetic_diversity/plot_Hs_saturation_curve_sampling_size_subgenlight_POPBYPOP.png"),
       width = 8, height = 5)




# ---- Method 2 : pops together, Hs ~ sampling size ----
## ---- downsample genlight ----
genlight_list <- list()
for (N in 1:30){
  samples_to_keep <- c()
  for (pop in levels(genlight$pop)) {
    # get samples name for the pop considered
    samples_all <- 
      data_samples %>% 
      dplyr::filter(id %in% genlight$ind.names) %>% 
      dplyr::filter(.data[[level]] == pop) %>% 
      dplyr::pull(id)
    
    
    # if more than N individuals in a pop, randomly sample N individuals names
    if(length(samples_all) > N) { 
      samples_temp <- sample(samples_all, N, replace=F) 
    }
    
    # if equal to N individuals in a pop, collect all individual names
    if(length(samples_all) <= N) {
      samples_temp <- samples_all 
    }
    
    samples_to_keep <- c(samples_to_keep, samples_temp)

  }
  
  genlight_list[[N]] <- 
    dartR::gl.keep.ind(x = genlight, 
                       ind.list = samples_to_keep,
                       mono.rm = TRUE)
}

genlight_list %>% 
  saveRDS("intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered_missind_callrate0.70_maf0.05_different_sampling_size_ALLPOPS.RDS")



## ---- iteratively compute Hs ----

het_list <- list()

for (N in 1:length(genlight_list)){
  het_list[[N]] <- 
    dartR::gl.report.heterozygosity(genlight_list[[N]], plot.out = F) %>% 
    mutate(N = N)
}

het_list %>% 
  saveRDS("intermediate/1_genetic_diversity/Alpha_diversity_dartR_different_sampling_size_ALLPOPS.RDS")

het_df <-
  het_list %>% 
  reduce(rbind)

het_df$site <- factor(x = het_df$pop, 
                      levels = levels(genlight@pop), 
                      ordered = TRUE)



## ---- plot saturation curve ----
ggplot(het_df, aes(nInd, He, color=site)) +
  scale_color_manual(values = color_perso) +
  geom_line(size = 0.8, alpha = 0.5) +
  xlim(0, 30)
ggsave(paste0("results/1_genetic_diversity/plot_Hs_saturation_curve_sampling_size_ALLPOPS.png"),
       width = 8, height = 5)




# 
# # ---- Resampling 100 times to N individuals ----
# 
# library(parallel)
# 
# N=6 # number of individuals by location
# i=1 # iteration
# 
# genlight_list <- lapply(1:100, function(i){
#   
#   cat(i, "/100 \n")
#   
#   samples_to_keep <- c()
#   for (pop in levels(genlight$pop)) {
#     
#     # get sample list from pop
#     samples_all <- 
#       data_samples %>% 
#       dplyr::filter(id %in% genlight$ind.names) %>% 
#       dplyr::filter(.[[level]] == pop) %>% 
#       dplyr::pull(id)
#     
#     # if more than N individuals in a population, sample N individuals names
#     if(length(samples_all) > N) { 
#       samples_sub <- sample(samples_all, N, replace=F) 
#     }
#     # if less than N individuals in a population, collect all individual names
#     if(length(samples_all) <= N) {
#       samples_sub <- samples_all 
#     }
#     
#     # add to sample list to keep
#     samples_to_keep <- c(samples_to_keep, samples_sub)
#   }
#   
#   # resample genlight
#   genlight_resampled <- 
#     dartR::gl.keep.ind(x = genlight, 
#                 ind.list = samples_to_keep,
#                 mono.rm = TRUE)
#   
#   # iteration
#   i = i+1
# })
# 
# genlight_list %>% 
#   saveRDS("intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered_missind_callrate0.70_maf0.05_resampled.RDS")

