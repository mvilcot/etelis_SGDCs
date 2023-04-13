## ---- load ----
dist_merge <-
  read_csv("results/3_distance_metrics/dist_geo_envt_res17-4_gd_sd.csv")

dist_mat <-
  readRDS("intermediate/3_distance_metrics/dist_geo_envt_res17-4_gd_sd.RDS")


## ---- parameters ----
# communities delineation
comm_delin = "taxonomic_scale_datasp2"
list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

names(dist_mat)

# parameters
comm = names(list_communities)[1]
metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"
metricDIST = "leastcost"


## ---- quick correlations ----
model <- cor(dist_merge[-c(1:3)])
corrplot::corrplot(model)


# ## ---- subset sites ----
# loc = "Seychelles"
# 
# for(i in 1:length(dist_mat)){
#   dist_mat[[i]] <- as.matrix(dist_mat[[i]])
#   dist_mat[[i]] <- dist_mat[[i]][!grepl(loc, rownames(dist_mat[[i]])),
#                                  !grepl(loc, colnames(dist_mat[[i]]))]
#   dist_mat[[i]] <- as.dist(dist_mat[[i]])
# }
# 
# dist_merge <- dist_merge[!grepl(loc, dist_merge$site),]



## ---- IBD simple ----

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

# # multiple
# MRM(dist_merge[[metricGD]] ~ 
#       dist_merge[[paste(names(list_communities)[1], "beta.jtu", sep = ".")]] + 
#       dist_merge[[paste(names(list_communities)[2], "beta.jtu", sep = ".")]] + 
#       dist_merge[[paste(names(list_communities)[3], "beta.jtu", sep = ".")]] , 
#     nperm = 9999)


## ---- plot ----
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
ggsave(width = 20, height = 6, 
       filename = paste0("results/4_continuity/IBD_SGDC_beta_allsites_", metricSD, "_", metricGD, "_", metricDIST, ".png"))




# ## ---- decomposition SGDCs ----
# 
# source("scripts/unsorted/_continuity/sgdcs_decomposition_Lamy.R")
# 
# SGDC.decomp(SD = dist_merge[[metricSD]], 
#             GD = dist_merge[[metricGD]], 
#             FACTOR = dist_merge[,c("environment","leastcost","seadist")])
# 




