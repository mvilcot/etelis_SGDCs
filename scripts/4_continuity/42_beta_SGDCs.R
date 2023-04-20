
## ---- parameters ----
# communities delineation
comm_delin = "taxonomic_scale_datasp2"
# comm_delin = "taxonomic_scale_Fishbase"
# comm_delin = "depth_category"
# comm_delin = "phylogenetic_distance"

list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

# parameters
comm = names(list_communities)[1]
metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"
metricDIST = "leastcost"


## ---- load ----
dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".csv"))

dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".RDS"))



## ---- subset sites ----
loc = "Seychelles"

for(i in 1:length(dist_mat)){
  dist_mat[[i]] <- as.matrix(dist_mat[[i]])
  dist_mat[[i]] <- dist_mat[[i]][!grepl(loc, rownames(dist_mat[[i]])),
                                 !grepl(loc, colnames(dist_mat[[i]]))]
  dist_mat[[i]] <- as.dist(dist_mat[[i]])
}

dist_merge <- dist_merge[!grepl(loc, dist_merge$site),]



## ---- 0. check assumptions ----
model <- cor(dist_merge[-c(1:3)])
corrplot::corrplot(model)

hist(dist_mat[[metricSD]])
hist(dist_mat[[metricGD]]) # GD not normal

shapiro.test(dist_mat[[metricSD]])
shapiro.test(dist_mat[[metricGD]]) # GD not normal



## ---- 1. IBD + SGDCs ----

# MRM
sMRM_IBDsd <- MRM(dist_mat[[metricSD]] ~ dist_mat[[metricDIST]], nperm = 9999)
sMRM_IBDgd <- MRM(dist_mat[[metricGD]] ~ dist_mat[[metricDIST]], nperm = 9999)
sMRM_SGDC <- MRM(dist_mat[[metricSD]] ~ dist_mat[[metricGD]], nperm = 9999)
# MRM(dist_merge[[metricSD]] ~ dist_merge[[metricDIST]], nperm = 9999) # same with distance matrix on long df

# Mantel
sMantel_IBDsd <- vegan::mantel(dist_mat[[metricSD]], dist_mat[[metricDIST]], permutations = 9999)
sMantel_IBDgd <- vegan::mantel(dist_mat[[metricGD]], dist_mat[[metricDIST]], permutations = 9999)
sMantel_SGDC <- vegan::mantel(dist_mat[[metricSD]], dist_mat[[metricGD]], permutations = 9999)

# distance decay models
sDecay_IBDsd <- decay.model(dist_mat[[metricSD]], dist_mat[[metricDIST]], y.type="dissim", model.type="pow", perm=999)
sDecay_IBDgd <- decay.model(dist_mat[[metricGD]], dist_mat[[metricDIST]], y.type="dissim", model.type="pow", perm=999)

plot(dist_mat[[metricDIST]], dist_mat[[metricSD]], pch = 16)
plot.decay(sDecay_IBDsd, col="green", remove.dots=T, add=T, lwd=2)

plot(dist_mat[[metricDIST]], dist_mat[[metricGD]], pch = 16)
plot.decay(sDecay_IBDgd, col="green", remove.dots=T, add=T, lwd=2)



## ---- plot ---- #
# add Seychelles color
dist_merge$Seychelles <- "No"
dist_merge[grep("Seychelles", dist_merge$site),]$Seychelles <- "Yes"

# species IBD
ggSD <- 
  ggplot(dist_merge) +
  geom_point(aes(.data[[metricDIST]], .data[[metricSD]], color = Seychelles), show.legend = F) +
  scale_color_brewer(palette="Dark2") +
  annotate('text', 
           x=min(dist_merge[[metricDIST]]), y=max(dist_merge[[metricSD]]),
           hjust = 0, vjust = 1,
           label=paste0("r Mantel = ", round(sMantel_IBDsd$statistic, 4), ", p = ", round(sMantel_IBDsd$signif, 5),
                        "\nr MRM = ", round(sMRM_IBDsd$r.squared[["R2"]], 4), ", p = ", round(sMRM_IBDsd$r.squared[["pval"]], 5))) +
  labs(title = "Species IBD")

# genetic IBD
ggGD <- 
  ggplot(dist_merge) +
  geom_point(aes(.data[[metricDIST]], .data[[metricGD]], color = Seychelles), show.legend = F) +
  scale_color_brewer(palette="Dark2") +
  annotate('text', 
           x=min(dist_merge[[metricDIST]]), y=max(dist_merge[[metricGD]]),
           hjust = 0, vjust = 1,
           label=paste0("r Mantel = ", round(sMantel_IBDgd$statistic, 4), ", p = ", round(sMantel_IBDgd$signif, 5),
                        "\nr MRM = ", round(sMRM_IBDgd$r.squared[["R2"]], 4), ", p = ", round(sMRM_IBDgd$r.squared[["pval"]], 5))) +
  labs(title = "Genetics IBD")


# SGDCs
ggSGDCs <-
  ggplot(dist_merge) +
  geom_point(aes(.data[[metricSD]], .data[[metricGD]], color = Seychelles)) +
  scale_color_brewer(palette="Dark2") +
  annotate('text', 
           x=min(dist_merge[[metricSD]]), y=max(dist_merge[[metricGD]]),
           hjust = 0, vjust = 1,
           label=paste0("r Mantel = ", round(sMantel_SGDC$statistic, 4), ", p = ", round(sMantel_SGDC$signif, 4),
                        "\nr MRM = ", round(sMRM_SGDC$r.squared[["R2"]], 4), ", p = ", round(sMRM_SGDC$r.squared[["pval"]], 5))) +
  labs(title = "SGDC")


# merge plots
ggSD + ggGD + ggSGDCs + plot_annotation(title = comm)
# ggsave(width = 20, height = 6, 
#        filename = paste0("results/4_continuity/beta_SGDC_IBD_allsites_", metricSD, "_", metricGD, "_", metricDIST, ".png"))





## ---- 2. SGDCs by community ----

gg_list <- list()
stat_list <- list()

comm = names(list_communities)[1]
for (comm in names(list_communities)){
  metricSDcomm = paste(comm, "beta.jtu", sep = ".")
  print(metricSDcomm)
  
  
  ## -- statistics --
  # Mantel test
  stat_list$Mantel[[comm]] <- vegan::mantel(dist_mat[[metricSDcomm]], dist_mat[[metricGD]], permutations = 9999)

  # MRM
  stat_list$MRM[[comm]] <- MRM(dist_mat[[metricSDcomm]] ~ dist_mat[[metricGD]], nperm = 9999)

  # LM
  stat_list$lm[[comm]] <- lm(dist_merge[[metricSDcomm]] ~ dist_merge[[metricGD]])

  
  ## -- plot --
  gg_list[[comm]] <- 
    ggplot(dist_merge, aes(.data[[metricSDcomm]], .data[[metricGD]])) +
    geom_point() +
    xlab(metricSDcomm) +
    # geom_text_repel(aes(label = sites), size=3) +
    # ggtitle(paste0(patt)) +
    annotate('text', 
             x=min(dist_merge[[metricSDcomm]]), y=max(dist_merge[[metricGD]]),
             hjust = 0, vjust = 1,
             label=paste0("r Mantel = ", round(stat_list$Mantel[[comm]]$statistic, 4), ", p = ", round(stat_list$Mantel[[comm]]$signif, 4),
                          "\nr MRM = ", round(stat_list$MRM[[comm]]$r.squared[["R2"]], 4), ", p = ", round(stat_list$MRM[[comm]]$r.squared[["pval"]], 5))) +
  labs(title = comm)
    
  
  # ggsave(gg, device = "png", height = 4, width = 6, units = "in",
  #        filename = paste0("results/04_continuity/SGDCs_beta_", patt, "_", level, "_", metricGD, "_", metricSD, "_", comm, ".png"))
  
}

gg_grob <- arrangeGrob(grobs = gg_list, ncol=4)
plot(gg_grob)
ggsave(gg_grob, width = 20, height = 10, 
       filename = paste0("results/4_continuity/beta_SGDCs_noSeychelles_", level, "_", metricGD, "_", metricSD, "_",  comm_delin, ".png"))




## ---- 3. MRM decomposition ----

names(dist_merge)
comm = names(list_communities)[1]

metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"

MRM(dist_mat[[metricGD]] ~ dist_mat$environment + dist_mat$seadist, nperm = 9999)
MRM(dist_mat[[metricSD]] ~ dist_mat$environment + dist_mat$seadist, nperm = 9999)

model <- cor(dist_merge[c("geodist",
               "seadist",
               "leastcost",
               "environment")])
corrplot::corrplot(model)




## ---- 4. decomposition variance ----

source("scripts/4_continuity/sgdcs_decomposition_Lamy.R")

comm = names(list_communities)[1]
metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"


png("results/4_continuity/variance_decomposition_beta_SGDCs.png",
    width = 10, height = 5, units = 'in', res = 300)
par(oma = c(2,12,2,2),
    mar = c(1,1,1,1))
SGDC.decomp(SD = dist_merge[[metricSD]],
            GD = dist_merge[[metricGD]],
            FACTOR = dist_merge[,c("environment","seadist")])
dev.off()






## ---- 5. SGDCs by distance phylo ----

SGDCs <- data.frame(community = names(list_communities))


temp <- lapply(SGDCs$community, FUN = function(x){stat_list$Mantel[[x]]$statistic})
names(temp) <- names(list_communities)

SGDCs <- 
  temp %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "community") %>% 
  dplyr::rename(rMantel = V1)

SGDCs$phylodist <- 
  gsub(pattern = "phylodist_", replacement = "", SGDCs$community)

plot(SGDCs$phylodist, SGDCs$rMantel)
summary(lm(SGDCs$phylodist ~ SGDCs$rMantel))




