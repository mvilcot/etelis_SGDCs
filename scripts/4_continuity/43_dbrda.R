# ---- read distance matrix ----
filters = "missind1_callrate0.70_maf0.05"
level = "site"
sites = "noSeychelles"
comm_delin = "taxonomy"
metric <- "Fst" # Lutjanidae.beta.jtu

# genlight <- 
#   read.genlight(filters, level,
#                 site2drop = "Seychelles", # NULL
#                 site2keep = NULL,
#                 station2drop = NULL,
#                 station2keep = NULL)

dist_merge <- read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))
dist_mat <- readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))

envt_site <- 
  read_csv("intermediate/3_distance_metrics/bio_oracle_variables_bdmean_stationbuffer_siteaverage.csv") # %>%
  # left_join(data_sites, by = "site")


# ---- subset Seychelles ----
if (sites == "noSeychelles"){
  envt_site <- mat.subset(envt_site, "Seychelles")
  dist_merge <- mat.subset(dist_merge, "Seychelles")
  dist_mat <- lapply(dist_mat, function(x) mat.subset(x, "Seychelles"))
  data_sites <- mat.subset(data_sites, "Seychelles")
}



# ---- 0. GENOMIC DATA (Y) ----
# ## change format ----
# genind <- dartR::gl2gi(genlight)
# levels(genind@pop)
# 
# 
# ## euclidean distances ----
# distgenEUCL <- dist(genind, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
# distgenEUCL <- dist_merge[[metric]]




# ---- 1. LON-LAT MEM (X) ----
## dbMEM ----
# https://github.com/laurabenestan/Seascape_reservebenefit/blob/main/04-db-rda/script-dbRDA-neutral.R

## in-water dbMEM ----
coord_sites <-
  data_sites %>%
  column_to_rownames("site") %>%
  dplyr::select(-number_samples)

# compute MEM
dbMEM_inwater <- adespatial::dbmem(dist_mat$seadist, MEM.autocor = "non-null", store.listw = TRUE)
dbMEM.vectors.inwater <- as.data.frame(dbMEM_inwater)

# Look at general output
summary(dbMEM_inwater)
dbmem <- as.data.frame(dbMEM_inwater)

## euclidian dist MEM ----

# library(adespatial)
# library(SoDA)
# 

# 
# # spatial coordinates to cartesian coordinates
# coord_sites_cart <- SoDA::geoXY(data_sites$latitude, data_sites$longitude, unit=1000)
# 
# # euclidean distance matrix
# euclidean_distances <- dist(coord_sites_cart, method="euclidean") 
# 
# # Perform dbMEM on euclidean distances
# dbMEM_eucl <- adespatial::dbmem(euclidean_distances, MEM.autocor = "non-null", store.listw = TRUE)
# 
# # Keeping default setting for truncation (length of the longest edge of the minimum spanning tree will be used as the threshold) and just positive MEM
# dbMEM_eucl = adespatial::dbmem(coord_sites) 
# 
# # Transform into a dataframe
# dbMEM.vectors.euclidian <- as.data.frame(dbMEM_eucl)


## choose between dbMEM inwater and dbMEM from euclidian distances ----

# cor.test(dbMEM.vectors.euclidian$MEM1, dbMEM.vectors.inwater$MEM1)







# ---- 2. ENVIRONMENT (X) ----
library(factoextra)

## scale ----
envt_site_scaled <-
  envt_site %>% 
  column_to_rownames("site") %>% 
  scale(center = T, scale = T) %>% 
  as.data.frame()

## PCA ----
envt_PCA <- prcomp(envt_site_scaled)
envt_PCAaxis <- as.data.frame(envt_PCA$x) # The principal components can be found in the $x matrix
fviz_eig(envt_PCA)

## subset the variables Env file ----
# pca_axis <- pca_axis[588:1055,]






# ---- 3. Run RDA ----
## all variables ----
# correlations between X axes
# X = cbind(dbMEM_inwater[,1:3], envt_PCAaxis[,1:4])
X = cbind(dbMEM_inwater[,1:2], envt_site_scaled[c(1,2,3,4,6)])

# genetic IBD/IBE
dbrda_GD0 <- vegan::capscale(dist_mat$Fst ~ 1, data = X) # null model, only intercept
dbrda_GDgeo <- vegan::capscale(dist_mat$Fst ~ ., data = dbMEM_inwater)
dbrda_GDenv <- vegan::capscale(dist_mat$Fst ~ ., data = envt_PCAaxis)
dbrda_GDall <- vegan::capscale(dist_mat$Fst ~ .,
                               data = X) # all variables

# species IBD/IBE
dbrda_SD0 <- vegan::capscale(dist_mat$Lutjanidae.beta.jtu ~ 1, data = X) # null model, only intercept
dbrda_SDgeo <- vegan::capscale(dist_mat$Lutjanidae.beta.jtu ~ MEM1 + MEM2 + MEM3 + MEM4 + MEM5, data = dbMEM_inwater)
dbrda_SDenv <- vegan::capscale(dist_mat$Lutjanidae.beta.jtu ~ PC1 + PC2 + PC3 + PC4 + PC5, data = envt_PCAaxis)
dbrda_SDall <- vegan::capscale(dist_mat$Lutjanidae.beta.jtu ~ ., 
                               data = X) # all variables

# SGDC
pcoaGD <- dudi.pco(quasieuclid(dist_mat$Fst), scannf=FALSE, nf=8)
dbrda_SGDC0 <- vegan::capscale(dist_mat$Lutjanidae.beta.jtu ~ 1, data = pcoaGD$li) 
dbrda_SGDC <- vegan::capscale(dist_mat$Lutjanidae.beta.jtu ~ ., data = pcoaGD$li) 

pcoaSD <- dudi.pco(quasieuclid(dist_mat$Lutjanidae.beta.jtu), scannf=FALSE, nf=8)
dbrda_GSDC0 <- vegan::capscale(dist_mat$Fst ~ 1, data = pcoaSD$li) 
dbrda_GSDC <- vegan::capscale(dist_mat$Fst ~ ., data = pcoaSD$li) 


## selection of explanatory variables ----

### anova ----
selGD <- ordiR2step(dbrda_GD0, scope = formula(dbrda_GDall), direction="both") 
selGD$anova

selSD <- ordiR2step(dbrda_SD0, scope = formula(dbrda_SDall), direction="both") 
selSD$anova

selSGDC <- ordiR2step(dbrda_SGDC0, scope = formula(dbrda_SGDC), direction="both") 
selSGDC$anova

selGSDC <- ordiR2step(dbrda_GSDC0, scope = formula(dbrda_GSDC), direction="both") 
selGSDC$anova

# ade4::s.value(dfxy = coord_sites, z = X[,c("MEM1")])
# ade4::s.value(dfxy = coord_sites, z = X[,c("MEM3")])


### correlation ----
library(car)
corX <- cor(X)
corrplot::corrplot(corX)

corENV <- cor(envt_site_scaled)
corrplot::corrplot(corENV)

vif(dbrda_SDall)




## final models ----
dbrda_GDfin <- vegan::capscale(dist_mat$Fst ~ MEM1 + BO2_dissoxmean_bdmean, data = X) # MEM 3
RsquareAdj(dbrda_GDfin)
anova(dbrda_GDfin)    

dbrda_SDfin <- vegan::capscale(dist_mat$Lutjanidae.beta.jtu ~ MEM1, data = X) 
RsquareAdj(dbrda_SDfin)
anova(dbrda_SDfin)    

dbrda_SGDCfin <- vegan::capscale(dist_mat$Lutjanidae.beta.jtu ~ A1, data = pcoaGD$li) 
RsquareAdj(dbrda_SGDCfin)
anova(dbrda_SGDCfin)    

dbrda_GSDCfin <- vegan::capscale(dist_mat$Fst ~ A1 + A2, data = pcoaSD$li) 
RsquareAdj(dbrda_SGDCfin) # same, no matter GD ~ SD or SD ~ GD
anova(dbrda_SGDCfin)    





# ---- 4. Plot ----
## GD ----
plot(dbrda_GDfin)

# get data
dbrda_GDsites <- 
  as.data.frame(vegan::scores(dbrda_GDall)$sites) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("site") %>% 
  left_join(data_sites, by = "site")

dbrda_GDvar <- 
  as.data.frame(dbrda_GDall$CCA$biplot[, 1:2]) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("variable")
  
ggplot() +
  # geom_encircle(aes(group = Protection,linetype = Protection,fill= Protection), s_shape = 1, expand = 0,
  #               alpha = 0.4, show.legend = FALSE) + # hull area 
  geom_point(data = dbrda_GDsites, aes(x= CAP1, y = CAP2, color = site)) +
  geom_text_repel(data = dbrda_GDsites, aes(x= CAP1, y = CAP2, color = site),
                  label = dbrda_GDsites$site, bg.color = "grey70", bg.r = 0.02) +
  scale_color_manual(values = color_perso) +
  geom_segment(data= dbrda_GDvar, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc")), show.legend = F) +
  geom_label_repel(data= dbrda_GDvar, 
                   aes(x= CAP1, y=CAP2, fontface=3),
                   label = dbrda_GDvar$variable,
                   label.size = NA,
                   size = 4,
                   fill = NA) +
  ggtitle(metric) +
  theme_light() +
  theme(legend.position="none")



## SD ----
plot(dbrda_SDfin)

# get data
dbrda_SDsites <- 
  as.data.frame(vegan::scores(dbrda_SDall)$sites) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("site") %>% 
  left_join(data_sites, by = "site")

dbrda_SDvar <- 
  as.data.frame(dbrda_SDall$CCA$biplot[, 1:2]) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("variable")

ggplot() +
  # geom_encircle(aes(group = Protection,linetype = Protection,fill= Protection), s_shape = 1, expand = 0,
  #               alpha = 0.4, show.legend = FALSE) + # hull area 
  geom_point(data = dbrda_SDsites, aes(x= CAP1, y = CAP2, color = site)) +
  geom_text_repel(data = dbrda_SDsites, aes(x= CAP1, y = CAP2, color = site),
                  label = dbrda_SDsites$site, bg.color = "grey70", bg.r = 0.02) +
  scale_color_manual(values = color_perso) +
  geom_segment(data= dbrda_SDvar, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc")), show.legend = F) +
  geom_label_repel(data= dbrda_SDvar, 
                   aes(x= CAP1, y=CAP2, fontface=3),
                   label = dbrda_SDvar$variable,
                   label.size = NA,
                   size = 4,
                   fill = NA) +
  ggtitle(metric) +
  theme_light() +
  theme(legend.position="none")




## SGDC ----
plot(dbrda_SGDCfin)
plot(dbrda_GSDCfin)

# # get data
# dbrda_SGDCsites <- 
#   as.data.frame(vegan::scores(dbrda_SGDCfin)$sites) %>% 
#   as_tibble(rownames = NA) %>% 
#   rownames_to_column("site") %>% 
#   left_join(data_sites, by = "site")
# 
# dbrda_SGDCvar <- 
#   as.data.frame(dbrda_SGDCfin$CCA$biplot) %>% 
#   as_tibble(rownames = NA) %>% 
#   rownames_to_column("variable")
# 
# ggplot() +
#   # geom_encircle(aes(group = Protection,linetype = Protection,fill= Protection), s_shape = 1, expand = 0,
#   #               alpha = 0.4, show.legend = FALSE) + # hull area 
#   geom_point(data = dbrda_SGDCsites, aes(x= CAP1, y = MDS1, color = site)) +
#   geom_text_repel(data = dbrda_SGDCsites, aes(x= CAP1, y = MDS1, color = site),
#                   label = dbrda_SGDCsites$site, bg.color = "grey70", bg.r = 0.02) +
#   scale_color_manual(values = color_perso) +
#   geom_segment(data= dbrda_SGDCvar, aes(x=0, xend=CAP1, y = 0, yend=CAP2), col = "black",
#                arrow=arrow(length=unit(0.01,"npc")), show.legend = F) +
#   geom_label_repel(data= dbrda_SGDCvar, 
#                    aes(x= CAP1, y=CAP2, fontface=3),
#                    label = dbrda_SGDCvar$variable,
#                    label.size = NA,
#                    size = 4,
#                    fill = NA) +
#   ggtitle(metric) +
#   theme_light() +
#   theme(legend.position="none")



# 5. Var decomp ----
vpGD <- varpart(dist_mat$Fst, dbMEM_inwater[,1:2], envt_PCAaxis[,1:2])
vpGD
plot(vpGD, bg = c(3, 5), Xnames = c("spatial", "environment"))

vpSD <- varpart(dist_mat$Lutjanidae.beta.jtu, dbMEM_inwater[,1:2], envt_PCAaxis[,1:3])
vpSD
plot(vpSD, bg = c(3, 5), Xnames = c("spatial", "environment"))



