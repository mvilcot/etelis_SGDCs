## ---- parameters ----
# communities delineation
comm_delin = "taxonomic_scale_datasp2"
list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

# parameters
comm = names(list_communities)[1]
metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"


## ---- quick correlations ----
model <- cor(dist_merge[-c(1:2)])
corrplot::corrplot(model)



## ---- statistics ----

# metricDIST = 'geodist'; matDIST = mat_geodist
metricDIST = 'lcdist'; matDIST = mat_lcdist

# LM
summary(lm(merge_beta[[metricSD]] ~ merge_beta[[metricDIST]]))
summary(lm(merge_beta[[metricGD]] ~ merge_beta[[metricDIST]]))

# MRM
MRM(merge_beta[[metricSD]] ~ merge_beta[[metricDIST]], nperm = 9999)
MRM(merge_beta[[metricGD]] ~ merge_beta[[metricDIST]], nperm = 9999)

# Mantel
stat_SDmantel <- vegan::mantel(as.dist(mat_SDbeta), as.dist(matDIST), permutations = 9999)
stat_GDmantel <- vegan::mantel(as.dist(mat_GDbeta), as.dist(matDIST), permutations = 9999)
stat_SGDCmantel <- vegan::mantel(as.dist(mat_SDbeta), as.dist(mat_GDbeta), permutations = 9999)



## ---- plot ----
# add Seychelles color
merge_beta$Seychelles <- "No"
merge_beta[grep("Seychelles", merge_beta[[level]]),]$Seychelles <- "Yes"

# species IBD
ggSD <- 
  ggplot(merge_beta) +
  geom_point(aes(.data[[metricDIST]], .data[[metricSD]], 
                 color = Seychelles), show.legend = F) +
  scale_color_brewer(palette="Dark2") +
  annotate('text', 
           x=min(merge_beta[[metricDIST]]), y=max(merge_beta[[metricSD]]),
           hjust = 0, vjust = 1,
           label=paste0("r Mantel = ", round(stat_SDmantel$statistic, 4), "\np = ", round(stat_SDmantel$signif, 5))) + 
  labs(title = "Species IBD")

# genetic IBD
ggGD <- 
  ggplot(merge_beta) +
  geom_point(aes(.data[[metricDIST]], .data[[metricGD]], 
                 color = Seychelles), show.legend = F) +
  scale_color_brewer(palette="Dark2") +
  annotate('text', 
           x=min(merge_beta[[metricDIST]]), y=max(merge_beta[[metricGD]]),
           hjust = 0, vjust = 1,
           label=paste0("r Mantel = ", round(stat_GDmantel$statistic, 4), "\np = ", round(stat_GDmantel$signif, 5))) +
  labs(title = "Genetics IBD")


# SGDCs
ggSGDCs <-
  ggplot(merge_beta) +
  geom_point(aes(.data[[metricSD]], .data[[metricGD]], color = Seychelles)) +
  scale_color_brewer(palette="Dark2") +
  # xlab(paste0(metricSD, " (", comm, ")")) +
  # geom_text_repel(aes(label = sites), size=3) +
  # ggtitle(paste0(patt)) +
  annotate('text', 
           x=min(merge_beta[[metricSD]]), y=max(merge_beta[[metricGD]]),
           hjust = 0, vjust = 1,
           label=paste0("r Mantel = ", round(stat_SGDCmantel$statistic, 4), "\np = ", round(stat_SGDCmantel$signif, 4))) +
  labs(title = "SGDC")


# merge plots
ggSD + ggGD + ggSGDCs + plot_annotation(title = comm)
ggsave(width = 20, height = 6, 
       filename = paste0("results/3_distance_decay/IBD_beta_noSeychelles_", level, "_", comm, "_", metricGD, "_", metricSD, "_", metricDIST, ".png"))




## ---- decomposition SGDCs ----

source("scripts/4_continuity/sgdcs_decomposition_Lamy.R")

SGDC.decomp(SD = merge_beta$beta.jtu, 
            GD = merge_beta$Fst, 
            FACTOR = merge_beta[,c("geodist","lcdist")])





