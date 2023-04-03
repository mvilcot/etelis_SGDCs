## ---- parameters ----
# filters
filters = "missind_callrate0.70_maf0.05"
level = "site"

# read genlight and remove locations with less than 2 individuals
genlight <- 
  read.genlight(filters, level, removeless2ind = TRUE)




## ---- Resampling 100 times to N individuals ----

library(parallel)

N=6 # number of individuals by location
i=1 # iteration

genlight_list <- lapply(1:100, function(i){
  
  cat(i, "/100 \n")
  
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
  
  # resample genlight
  genlight_resampled <- 
    gl.keep.ind(x = genlight, 
                ind.list = samples_to_keep,
                mono.rm = TRUE)
  
  i = i+1
})

saveRDS(genlight_list, 
        "intermediate/01_genetic_diversity/Genlight_Etelis_coruscans_ordered_missind_callrate0.70_maf0.05_resampled.RDS")

