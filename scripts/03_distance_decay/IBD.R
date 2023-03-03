
# library(adegenet)
# library(hierfstat)
# library(ecodist)
# library(vegan)
# library(reshape2)

# https://popgen.nescent.org/StartSNP.html
# https://popgen.nescent.org/DifferentiationSNP.html



## ---- load data ----

# read GD and SD genetic diversity
level = "site"

gd_beta <- readRDS(paste0("results/01_genetic_diversity/gd_list_pairwise_", level, ".RDS"))
sd_beta <- readRDS(paste0("results/02_species_diversity/sd_list_pairwise_", level, ".RDS"))
# sd_beta <- readRDS(paste0("results/02_species_diversity/sd_list_pairwise_", level, "_phylogenetic_scale.RDS"))


# communities delineation
list_communities <- readRDS("intermediate/02_species_diversity/List_community.RDS")
# list_communities <- readRDS("intermediate/02_species_diversity/List_community_phylogenetic_scale.RDS")

# parameters
metricSD = "beta.jtu"
metricGD = "Fst"
comm = names(list_communities)[1]


## ---- setup beta data ----
# get distance matrix
mat_SDbeta <- as.matrix(sd_beta[[comm]][[metricSD]])
mat_GDbeta <- as.matrix(gd_beta[[metricGD]])

# keep only sampling sites present in both GD and SD
mat_SDbeta <- mat_SDbeta[rownames(mat_SDbeta) %in% rownames(mat_GDbeta),
                         colnames(mat_SDbeta) %in% colnames(mat_GDbeta)]
mat_GDbeta <- mat_GDbeta[rownames(mat_GDbeta) %in% rownames(mat_SDbeta),
                         colnames(mat_GDbeta) %in% colnames(mat_SDbeta)]

# order rows alphabetically
mat_SDbeta <- mat_SDbeta[order(rownames(mat_SDbeta)), order(colnames(mat_SDbeta))]
mat_GDbeta <- mat_GDbeta[order(rownames(mat_GDbeta)), order(colnames(mat_GDbeta))]

# pivot longer distance matrix
melt_SDbeta <- melt.dist(dist = mat_SDbeta, metric = metricSD)
melt_GDbeta <- melt.dist(dist = mat_GDbeta, metric = metricGD)

# merge two distance matrix into one df
merge_beta <- 
  left_join(melt_SDbeta, melt_GDbeta,
            by = c(paste0(level, "1"), paste0(level, "2")))

# create one column for both locations info
merge_beta[[level]] <- paste0(merge_beta[[paste0(level, "1")]], 
                              "-", 
                              merge_beta[[paste0(level, "2")]])

# relocate to first column
merge_beta <- 
  merge_beta %>% relocate(level)




## ---- geographic distance ----

# get mean coordinates by location
data_coord <- 
  data_sites %>% 
  # dplyr::select(all_of(level), Longitude_approx, Latitude_approx) %>% 
  # dplyr::rename("longitude" = "Longitude_approx",
  #        "latitude" = "Latitude_approx") %>%
  group_by(.data[[level]]) %>% 
  summarise(longitude=mean(Longitude_approx),
          latitude=mean(Latitude_approx)) %>% 
  column_to_rownames(level)
  
# compute distance between locations
mat_geodist <- geodist(data_coord, measure = "geodesic")
colnames(mat_geodist) <- rownames(mat_geodist) <- rownames(data_coord)
mat_geodist <- mat_geodist[order(rownames(mat_geodist)), order(colnames(mat_geodist))]

# keep only sites present in SD and GD
mat_geodist <- mat_geodist[rownames(mat_geodist) %in% rownames(mat_SDbeta),
                         colnames(mat_geodist) %in% colnames(mat_SDbeta)]

# pivot longer distance matrix
melt_geodist <- melt.dist(dist = mat_geodist, metric = "geodist")

# merge two distance matrix into one df
merge_beta <- 
  merge_beta %>% 
  left_join(melt_geodist,
            by = c(paste0(level, "1"), paste0(level, "2"))) 




## ---- least-cost distance ----

# get bathymetry data
Bathy <- getNOAA.bathy(lon1 = -180, lon2 = 180,
                       lat1 = -25, lat2 = 25,
                       resolution = 100)
saveRDS(Bathy, "intermediate/03_distance_decay/bathymetry.RDS")
color_blues <- colorRampPalette(c("purple", "blue", "cadetblue1", "cadetblue2", "white"))

# save bathymetry map
pdf(file = "results/03_distance_decay/bathymetry_plot.pdf", height=7, width=14)
plot.bathy(Bathy, step=2000, deepest.isobath = -14000, shallowest.isobath = 0, image = TRUE, bpal = color_blues(20))
scaleBathy(Bathy, deg = 5, x = "topleft", inset = 5)
points(data_coord[["longitude"]], data_coord[["latitude"]],
       pch = 21, col = "black", bg = "red", cex = 1.3)
dev.off()

# compute least-cost distance
trans <- trans.mat(Bathy, min.depth=-5, max.depth=NULL)
mat_lcdist <- 
  lc.dist(trans, data_coord, res = "dist") %>% 
  as.matrix()
# write.csv(as.matrix(dist_geo), file = "intermediate/03_distance_decay/least_cost_distance.csv")

# keep only sites present in SD and GD
mat_lcdist <- mat_lcdist[order(rownames(mat_lcdist)), order(colnames(mat_lcdist))]
mat_lcdist <- mat_lcdist[rownames(mat_lcdist) %in% rownames(mat_SDbeta),
                           colnames(mat_lcdist) %in% colnames(mat_SDbeta)]

# pivot longer distance matrix
melt_lcdist <- melt.dist(dist = mat_lcdist, metric = "lcdist")

# merge two distance matrix into one df
merge_beta <- 
  merge_beta %>% 
  left_join(melt_lcdist,
            by = c(paste0(level, "1"), paste0(level, "2"))) 





## ---- subset to specific locations ----
patt = "Seychelles"

mat_SDbeta <- mat_SDbeta[!grepl(patt, rownames(mat_SDbeta)),
                         !grepl(patt, colnames(mat_SDbeta))]
mat_GDbeta <- mat_GDbeta[!grepl(patt, rownames(mat_GDbeta)),
                         !grepl(patt, colnames(mat_GDbeta))]
mat_geodist <- mat_geodist[!grepl(patt, rownames(mat_geodist)),
                           !grepl(patt, colnames(mat_geodist))]
mat_lcdist <- mat_lcdist[!grepl(patt, rownames(mat_lcdist)),
                         !grepl(patt, colnames(mat_lcdist))]

merge_beta <- merge_beta[!grepl(patt, merge_beta[[level]]),]






## ---- statistics ----

# metricDIST = 'geodist'; matDIST = mat_geodist
metricDIST = 'lcdist'; matDIST = mat_lcdist


# # LM
# summary(lm(merge_beta[[metricSD]] ~ merge_beta[[metricDIST]]))
# summary(lm(merge_beta[[metricGD]] ~ merge_beta[[metricDIST]]))
# 
# 
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
  ggplot(merge_beta, ) +
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
       filename = paste0("results/03_distance_decay/IBD_beta_noSeychelles_", level, "_", comm, "_", metricGD, "_", metricSD, "_", metricDIST, ".png"))


## ---- distance decay ----

mat_SDbeta <- as.dist(mat_SDbeta)
mat_GDbeta <- as.dist(mat_GDbeta)
mat_lcdist <- as.dist(mat_lcdist)

# distance decay models
decaySD_exp <- decay.model(mat_SDbeta, mat_lcdist, y.type="dissim", model.type="exp", perm=999)
decaySD_pow <- decay.model(mat_SDbeta, mat_lcdist, y.type="dissim", model.type="pow", perm=999)
decayGD_exp <- decay.model(mat_GDbeta, mat_lcdist, y.type="dissim", model.type="exp", perm=999)
decayGD_pow <- decay.model(mat_GDbeta, mat_lcdist, y.type="dissim", model.type="pow", perm=999)

# plots
png(paste0("results/03_distance_decay/Decay_beta_noSeychelles_", level, "_", comm, "_", metricGD, "_", metricSD, "_", metricDIST, ".png"),
    width = 15, height = 6, units = 'in', res = 300)
par(mfrow=c(1,2))

# plot(rnorm(100),typ="n",xlim=c(0,10000),ylim=c(0,0.5),cex.lab=0.6,cex.axis=0.5,mgp=c(2,0.5,0),tcl=-0.25,
#      xlab="",ylab="",las=1,axes=F)
# mtext(text = "\u03b2 diversity between sample sites",side=2,line=1.7,at=0.25,cex=0.7) 
# mtext(text ="Geographical distance between sample sites (km)",side=1,line=1.5,at=5000,cex=0.7)
# mtext(text="(b)",side=2,line=1.7,at=0.6,las=2,cex=1.2,font=3)
# 
# axis(side=1,at=c(0,2000,4000,6000,8000,10000),labels =c("0","2000","4000","6000","8000","10000"),
#      tcl=-0.25,cex.axis=0.6,mgp=c(3,0.25,0))
# axis(side=2,at=seq(0,0.5,0.1),labels=seq(0,0.5,0.1),tcl=-0.25,cex.axis=0.6,las=2,mgp=c(3,0.5,0))

plot(mat_lcdist, mat_SDbeta, pch = 16)
plot.decay(decaySD_exp, col="green", remove.dots=T, add=T, lwd=2)
plot.decay(decaySD_pow, col="blue", remove.dots=T, add=T, lwd=2)

plot(mat_lcdist, mat_GDbeta, pch = 16)
plot.decay(decayGD_exp, col="green", remove.dots=T, add=T, lwd=2)
plot.decay(decayGD_pow, col="blue", remove.dots=T, add=T, lwd=2)

dev.off()

# stats values

stats_decay <- list(decaySD_exp = decaySD_exp,
                    decaySD_pow = decaySD_pow,
                    decayGD_exp = decayGD_exp,
                    decayGD_pow = decayGD_pow)

table_decay <- do.call(rbind,lapply(1:length(stats_decay),function(i){
  c(names(stats_decay)[i],stats_decay[[i]]$b.slope,stats_decay[[i]]$pseudo.r.squared,stats_decay[[i]]$p.value)
}))
colnames(table_decay) <- c("model","Slope","Pseudo_R2","P.value")
table_decay <- data.frame(model=table_decay[,1],
                    apply(table_decay[,2:4],2,as.numeric))



