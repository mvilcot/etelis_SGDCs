# ---- Setup ----
## parameters ----
## communtity delineation
comm_delin_list <-
  c("taxonomy",
    "taxonomy_depth1_crosses45-400m",
    "taxonomy_depth2_contains45-400m",
    "taxonomy_depth3_within45-400m",
    "taxonomy_env_reef-associated")

comm_delin <- comm_delin_list[1]

# list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))
names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")

sites = "noSeychelles"


## load ----
dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))

envt_site <- 
  read_csv("intermediate/3_distance_metrics/bio_oracle_variables_bdmean_stationbuffer_siteaverage.csv") # %>%


## subset sites ----
if (sites == "noSeychelles"){
  dist_mat <- lapply(dist_mat, function(x) mat.subset(x, "Seychelles"))
  dist_merge <- mat.subset(dist_merge, "Seychelles")
}



# ---- 0. check assumptions ----
comm <- "Lutjanidae"
metricSD <- paste(comm, "beta.jtu", sep = ".")
metricGD <- "Fst"
metricDIST <- "seadist"

corALL <- cor(dist_merge[-c(1:3)])
corrplot::corrplot(corALL)

hist(dist_mat[[metricSD]])
hist(dist_mat[[metricGD]]) # GD not normal

shapiro.test(dist_mat[[metricSD]])
shapiro.test(dist_mat[[metricGD]]) # GD not normal
qqnorm(dist_mat[[metricSD]])
qqnorm(dist_mat[[metricGD]])



# ---- 1. IBD or IBE ----
## lon-lat dbMEM (X) ----
# https://github.com/laurabenestan/Seascape_reservebenefit/blob/main/04-db-rda/script-dbRDA-neutral.R
# compute MEM
dbMEM_seadist <- adespatial::dbmem(dist_mat$seadist, store.listw = TRUE) # MEM.autocor = "non-null",
dbMEM.vectors.inwater <- as.data.frame(dbMEM_seadist)

# select MEM
dbrda_GD0 <- vegan::capscale(dist_mat$Fst ~ 1, data = dbMEM_seadist) # null model, only intercept
dbrda_GDgeo <- vegan::capscale(dist_mat$Fst ~ ., data = dbMEM_seadist)
selGEO <- ordiR2step(dbrda_GD0, scope = formula(dbrda_GDgeo), direction="both") 
selGEO$anova



## environment PCA (X) ----
library(factoextra)

# scale
envt_site_scaled <-
  envt_site %>% 
  column_to_rownames("site") %>% 
  scale(center = T, scale = T) %>% 
  as.data.frame()

# PCA
PCA_envt <- prcomp(envt_site_scaled)
PCA_envtaxis <- as.data.frame(PCA_envt$x) # The principal components can be found in the $x matrix
fviz_eig(PCA_envt)


## plot ----
metricGD <- "Fst"
metricSDs <- paste(names_communities, "beta.jtu", sep = ".")
metricDIST <- 'seadist' #'environment'

gg_list <- list()
i=1

library(lmPerm)

for (metricDIV in c(metricGD, metricSDs)){
  # for (metricDIST in c("seadist")){ #"environment"
  
  for (locations in c("allsites", "noSeychelles")){
    
    ## get locations
    dist_mergeSUB <- dist_merge
    dist_matSUB <- dist_mat
    if(locations == "noSeychelles"){
      dist_mergeSUB <- mat.subset(dist_merge, "Seychelles")
      dist_matSUB <- lapply(dist_mat, function(x) mat.subset(x, "Seychelles"))
    }
    
    ## add Seychelles color
    dist_mergeSUB$Seychelles <- "No"
    dist_mergeSUB[grep("Seychelles", dist_mergeSUB$site),]$Seychelles <- "Yes"
    
    
    ## Mantel
    sMantel <- vegan::mantel(dist_matSUB[[metricDIV]], dist_matSUB[[metricDIST]], permutations = 9999)

    if(metricDIST == 'seadist') color_x <- "#3182BD"
    if(metricDIST == 'environment') color_x <- "#31A354"
    if(metricDIV %in% metricGD) color_y <- "#F69C73FF" # rocket
    if(metricDIV %in% metricSDs) color_y <- "#395D9CFF" # mako
    
    # proper beta and families
    metricLAB <- gsub('beta.', 'β', metricDIV)
    if(metricDIV %in% metricSDs){
        temp <- str_split_1(metricLAB, '[.]')
      metricLAB <- paste0(temp[2], ' (', temp[1], ')')
    }

    
    if(sMantel$signif < 0.05){
      gg_list[[i]] <-
        ggplot(dist_mergeSUB, aes(.data[[metricDIST]], .data[[metricDIV]])) + #, color = Seychelles
        geom_smooth(method = "lm", color = "grey50", se = F) +
        geom_point(aes(shape = Seychelles), show.legend = F, color = color_y) +
        scale_shape_manual(values = c(19, 1)) +
        annotate('text',
                 x=min(dist_mergeSUB[[metricDIST]]), y=max(dist_mergeSUB[[metricDIV]]),
                 hjust = 0, vjust = 1,
                 # label=paste0("adj R² = ", sprintf("%.3f", sdbRDA$statistic), "\np = ", sdbRDA$signif)) +
                 label=paste0("r Mantel = ", sprintf("%.3f", sMantel$statistic), "\np = ", sMantel$signif)) +
        theme_light() +
        ylab(metricLAB) +
        theme(
          axis.title.x = element_text(colour = color_x),
          # axis.title.y = element_text(colour = color_y)
        )
    }
    
    if(sMantel$signif > 0.05){
      gg_list[[i]] <-
        ggplot(dist_mergeSUB, aes(.data[[metricDIST]], .data[[metricDIV]], shape = Seychelles)) + #
        geom_point(aes(shape = Seychelles), show.legend = F, color = color_y) +
        scale_shape_manual(values = c(19, 1)) +
        annotate('text',
                 x=min(dist_mergeSUB[[metricDIST]]), y=max(dist_mergeSUB[[metricDIV]]),
                 hjust = 0, vjust = 1,
                 label=paste0("adj R² = ", sprintf("%.3f", sdbRDA$statistic), "\np = ", sdbRDA$signif)) +
                 # label=paste0("r Mantel = ", sprintf("%.3f", sMantel$statistic), "\np = ", sMantel$signif)) +
        theme_light() +
        ylab(metricLAB) +
        theme(
          axis.title.x = element_text(colour = color_x),
          # axis.title.y = element_text(colour = color_y)
        )
    }
    
    i = i+1
  }
}


### {FIGURES S8 and S9} ####

gg_grob <-
  patchwork::wrap_plots(gg_list, ncol = 2, tag_level = "new") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
plot(gg_grob)
ggsave(gg_grob, width = 10, height = 15, dpi = 500,
       filename = paste0("results/4_continuity/_S9_isolation_by_", metricDIST, "_", comm_delin, ".png"))






# ---- 2. SGDC table ----
## table ----
library(boot)

metricGD = "Fst"

# create empty table

# coeff_function <- function(formula, data, indices) {
#   d <- data[indices,] #allows boot to select sample
#   fit <- lm(formula, data=d) #fit regression model
#   return(summary(fit)$coefficients[2,1]) #return R-squared of model
# }
library(ecodist)        # MRM

SGDC_beta <- tibble()
i=1 # initialisation
for(comm_delin in comm_delin_list[1]) { #[1]

  for (locations in c("all sites", "without Seychelles")){

    dist_mat <- readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))
    dist_merge <- read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

    ## subset sites --
    if(locations == "without Seychelles"){
      dist_mat <- lapply(dist_mat, function(x){mat.subset(x, "Seychelles")})
      dist_merge <- mat.subset(dist_merge, "Seychelles")
    }

    for (comm in names_communities){

      for (metricSD in c("beta.jac", "beta.jtu", "beta.jne")){

        metricSDcomm <- paste(comm, metricSD, sep = ".")

        ### Mantel ----
        sMantel <- vegan::mantel(dist_mat[[metricSDcomm]], dist_mat[[metricGD]], permutations = 9999)
        
        ## Get bootstraps on R²
        fm <- as.formula(paste("scale(", metricSDcomm, ") ~ scale(", metricGD, ")", sep=""))
        sMantel2 <- summary(lm(fm, data = dist_merge))
        # boots <- boot(data=dist_merge, statistic=coeff_function, R=1000, formula=fm)
        
        fm <- as.formula(paste(metricSDcomm, " ~ ", metricGD, sep=""))
        sMantel3 <- ecodist::mantel(dist_mat[[metricSDcomm]] ~ dist_mat[[metricGD]], nperm = 9999, nboot = 9999, cboot = 0.95)        
        # boots_ci <- boot.ci(boots, type="bca")
        # plot(reps)

        
        ### Pearson ----
        sPearson <- cor.test(dist_merge[[metricSDcomm]], dist_merge[[metricGD]],
                             alternative = c("two.sided"), method = "pearson")
        sPearsonPERM <- jmuOutlier::perm.cor.test(dist_merge[[metricSDcomm]], dist_merge[[metricGD]], 
                                                  alternative = c("two.sided"), method = "pearson", 
                                                  num.sim = 9999)
        
        
        ### Procruste ----
        sProtest <- vegan::protest(dist_merge[[metricGD]], dist_merge[[metricSDcomm]])
        
        
        ### dbRDA ----
        # pcoaSD <- prcomp(quasieuclid(dist_mat[[metricSDcomm]]))
        # # # pcoaSD <- prcomp(dist_mat[[metricSDcomm]])
        # # # factoextra::fviz_eig(pcoaSD)
        # # ## variable selection
        # # dbrda_SGDC0 <- vegan::capscale(dist_mat[[metricGD]] ~ 1, data = as.data.frame(pcoaSD$x[, 1:7]))
        # # dbrda_SGDCall <- vegan::capscale(dist_mat[[metricGD]] ~ ., data = as.data.frame(pcoaSD$x[, 1:7]))
        # # selSGDC <- ordistep(dbrda_SGDC0,
        # #                       scope = formula(dbrda_SGDCall),
        # #                       direction="both")
        # # 
        # # RsquareAdj(selSGDC) # rsquared RDA
        # # anova.cca(selSGDC) # global signicativity of the RDA
        # 
        # dbrda_SGDCfin <- vegan::capscale(dist_mat[[metricGD]] ~ ., data = as.data.frame(pcoaSD$x[, 1:7]))
        
        

        ### Store results in table ----
        SGDC_beta[i, "community_delineation"] <- comm_delin
        SGDC_beta[i, "locations"] <- locations
        SGDC_beta[i, "taxonomic_scale"] <- comm
        SGDC_beta[i, "metricGD"] <- metricGD
        SGDC_beta[i, "metricSD"] <- metricSD
        SGDC_beta[i, "Mantel_r"] <- sMantel$statistic #sMantel3[1]
        SGDC_beta[i, "Mantel_rlwr"] <- sMantel2$coefficients[2,1] - 1.96*sMantel2$coefficients[2,2] # or sMantel3[5]  # or boots_ci$bca[4] 
        SGDC_beta[i, "Mantel_rupr"] <- sMantel2$coefficients[2,1] + 1.96*sMantel2$coefficients[2,2] # or sMantel3[6]  # orboots_ci$bca[5]
        SGDC_beta[i, "Mantel_pval"] <- sMantel$signif #sMantel3[2] ##Ha: r>0
        # SGDC_beta[i, "dbRDA_R2"] <- RsquareAdj(dbrda_SGDCfin)$r.squared #selSGDC
        # SGDC_beta[i, "dbRDA_adjR2"] <- RsquareAdj(dbrda_SGDCfin)$adj.r.squared
        # SGDC_beta[i, "dbRDA_pval"] <- anova.cca(dbrda_SGDCfin)[1,4]
        SGDC_beta[i, "Pearson_r"] <- sPearson$estimate
        SGDC_beta[i, "Pearson_pval"] <- sPearsonPERM$p.value
        SGDC_beta[i, "Procruste_r"] <- sqrt(1-sProtest$ss)
        SGDC_beta[i, "Procruste_pval"] <- sProtest$signif
        
        cat(comm_delin, locations, comm, metricSD, "\n")
        i=i+1
      }
    }
  }
}


SGDC_beta$Mantel_signif <- ifelse(SGDC_beta$Mantel_pval < 0.001, "***", 
                                  ifelse(SGDC_beta$Mantel_pval < 0.01, "**", 
                                         ifelse(SGDC_beta$Mantel_pval < 0.05, "*", "NS")))
# SGDC_beta$dbRDA_signif <- ifelse(SGDC_beta$dbRDA_pval < 0.001, "***",
#                                  ifelse(SGDC_beta$dbRDA_pval < 0.01, "**",
#                                         ifelse(SGDC_beta$dbRDA_pval < 0.05, "*", "NS")))
SGDC_beta$Pearson_signif <- ifelse(SGDC_beta$Pearson_pval < 0.001, "***",
                                 ifelse(SGDC_beta$Pearson_pval < 0.01, "**",
                                        ifelse(SGDC_beta$Pearson_pval < 0.05, "*", "NS")))
SGDC_beta$Procruste_signif <- ifelse(SGDC_beta$Procruste_pval < 0.001, "***", 
                                   ifelse(SGDC_beta$Procruste_pval < 0.01, "**", 
                                          ifelse(SGDC_beta$Procruste_pval < 0.05, "*", "NS")))

# plot(SGDC_beta$Mantel_r^2, SGDC_beta$Procruste_r^2)
# cor.test(SGDC_beta$Mantel_r^2, SGDC_beta$Procruste_r^2)
# 
# temp <- 
#   SGDC_beta %>% 
#   filter(metricSD != "beta.jne") %>% 
#   mutate(LABEL = paste(locations, taxonomic_scale, metricSD, sep = '-'))
# gg <-
#   ggplot(temp, aes(Mantel_pval, Procruste_pval, label = LABEL)) +
#   geom_point() 
# plotly::ggplotly(gg)
# cor.test(temp$Mantel_pval, temp$Procruste_pval)


SGDC_beta %>%
  write_csv("results/4_continuity/beta_SGDCs_table_CIlm_Pearson_Procruste_taxonomy.csv")

# temp <- 
#   SGDC_beta %>%
#   dplyr::filter(metricSD != "beta.jne") %>% 
#   mutate(Mantel_pval = paste(round(Mantel_pval, 3), Mantel_signif)) %>% 
#   mutate(Procruste_pval = paste(round(Procruste_pval, 3), Procruste_signif)) %>% 
#   mutate(Mantel_pval = gsub(" NS", "", Mantel_pval)) %>% 
#   mutate(Procruste_pval = gsub(" NS", "", Procruste_pval)) %>% 
#   dplyr::select(locations, metricGD, metricSD, taxonomic_scale, Mantel_r, Mantel_pval, Procruste_r, Procruste_pval) %>% 
#   dplyr::arrange(locations, metricSD)
# 
# temp 
# temp %>% 
#   write_csv("results/4_continuity/_beta_SGDCs_table_CIlm_Pearson_Procruste_taxonomy.csv")





## plot one community delin ----
SGDC_betaSUB <-
  SGDC_beta %>% 
  dplyr::filter(community_delineation == "taxonomy") %>% 
  dplyr::filter(locations == "without Seychelles") %>% 
  dplyr::filter(metricSD != "beta.jne")

SGDC_betaSUB$taxonomic_scale <- 
  factor(SGDC_betaSUB$taxonomic_scale,
         level = names_communities)

SGDC_betaSUB$metricSD <- 
  factor(SGDC_betaSUB$metricSD,
         level = c("beta.jac", "beta.jtu")) #"beta.jne", 

ggplot(SGDC_betaSUB, aes(
                       x=Pearson_r,
                       y=taxonomic_scale,
                       shape=Pearson_signif,
                       group=metricSD,
                       color=metricSD)) +
  geom_vline(xintercept = 0, linetype=2, color="grey20") +
  geom_point(aes(x = Pearson_r)) +
  # geom_pointrange(aes(x=Mantel_r,
  #                     xmin=Mantel_rlwr,
  #                     xmax=Mantel_rupr),
  #                   size = 0.3, linewidth = 0.3,
  #                   position=position_dodge(width = 0.4)) +
  scale_shape_manual(values=c(19, 19, 1), )+
  scale_color_manual(values = c("#60CEACFF", "#395D9CFF"))+ # "0B0405FF", 
  # scale_color_viridis(option="mako", discrete=T, end = 0.8) +
  xlim(c(-0.6, 0.9)) +
  ylab("") +
  xlab("r") +
  theme_light() +
  theme(legend.title=element_blank(), legend.position="top") +
  guides(shape = "none") +
  scale_y_discrete(limits=rev)

ggsave(paste0("results/4_continuity/beta_SGDCs_DWplot_taxonomy_noSeychelles_Pearson.png"),
       height = 4, width = 5, units = 'in', dpi = 300)



## plot with all comm delin ----
SGDC_beta$taxonomic_scale <- factor(SGDC_beta$taxonomic_scale,
                                    level = names_communities)

ggplot(SGDC_beta, aes(x=community_delineation,
                      y=Pearson_r,
                      shape=Pearson_signif,
                      # y=dbRDA_adjR2,
                      # shape=dbRDA_signif,
                      group=taxonomic_scale,
                      color=taxonomic_scale)) +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(size = 2.5, position=position_dodge(width = 0.4)) +
  scale_shape_manual(values=c(19, 19, 19, 1))+
  scale_color_manual(values = c("#92e101", "#78c556", "#5fa8aa", "#458cff"))+
  facet_grid(vars(metricSD), vars(locations), scale="fixed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(paste0("results/4_continuity/beta_SGDCs_DWplot_Pearson.png"),
       height = 8, width = 8, units = 'in', dpi = 300)



# ---- 3. SGDCs plot ----
SGDC_beta <- read_csv("results/4_continuity/beta_SGDCs_table_CIlm_dbRDA_Pearson_taxonomy.csv")

SGDC_betaSUB <- 
  SGDC_beta %>% 
  # dplyr::filter(signif != "NS") %>% 
  dplyr::filter(community_delineation == "taxonomy") %>% 
  dplyr::filter(metricSD == "beta.jtu") %>% 
  dplyr::filter(locations == "without Seychelles")

comm_delin <- "taxonomy"
gg_list <- list()
metricGD = "Fst"

for (i in 1:nrow(SGDC_betaSUB)){
  comm <- SGDC_betaSUB[i, "taxonomic_scale"] %>% pull()
  metricSD <- SGDC_betaSUB[i, "metricSD"] %>% pull()
  metricSDcomm <- paste0(comm, ".", metricSD)
  
  metricXLAB <- gsub('beta.', 'β', metricSD)

  dist_merge <-
    read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))
  
  if(SGDC_betaSUB[i, "locations"] == "without Seychelles"){
    dist_merge <- mat.subset(dist_merge, "Seychelles")
  }

  if(SGDC_betaSUB[i, "Pearson_signif"] != "NS"){
    gg_list[[i]] <-
      ggplot(dist_merge, aes(.data[[metricSDcomm]], .data[[metricGD]])) +
      geom_smooth(method = "lm", color = "grey50", se = FALSE) +
      geom_point(shape = 19, alpha = 0.8) +
      xlab(paste0(metricXLAB, " (", comm, ")")) +
      annotate('text',
               x = min(dist_merge[[metricSDcomm]]), y = 1.1 * max(dist_merge[[metricGD]]),
               hjust = 0, vjust = 1,
               label = paste0("r = ", round(SGDC_betaSUB[i, "Pearson_r"], 4), "\np = ", round(SGDC_betaSUB[i, "Pearson_pval"], 4))) +
      theme_light()
  }
  if(SGDC_betaSUB[i, "Pearson_signif"] == "NS"){
    gg_list[[i]] <-
      ggplot(dist_merge, aes(.data[[metricSDcomm]], .data[[metricGD]])) +
      geom_point(shape = 19, alpha = 0.8) +
      xlab(paste0(metricXLAB, " (", comm, ")")) +
      annotate('text',
               x = min(dist_merge[[metricSDcomm]]), y = 1.1 * max(dist_merge[[metricGD]]),
               hjust = 0, vjust = 1,
               label = paste0("r = ", round(SGDC_betaSUB[i, "Pearson_r"], 4), "\np = ", round(SGDC_betaSUB[i, "Pearson_pval"], 4))) +
      theme_light()
  }
  
}

### {FIGURE 4} ####
gg_grob <-
  patchwork::wrap_plots(gg_list, tag_level = "new") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
plot(gg_grob)

ggsave(gg_grob, width = 8, height = 7, dpi = 500,
       filename = paste0("results/4_continuity/_4_beta_SGDCs_plot_significant_REVISION_Pearson.png"))
ggsave(gg_grob, width = 8, height = 7, 
       filename = paste0("results/4_continuity/_4_beta_SGDCs_plot_significant_REVISION_Pearson.pdf"))







# ---- 4. Variance decomposition----

## load ----
source("scripts/4_continuity/sgdcs_decomposition_Lamy.R")
## >>>> to adapt to beta diversity... with phytools::multi.mantel? ####

comm = names_communities[2]
metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"

## subset Seychelles ----
dist_merge <- mat.subset(dist_merge, "Seychelles")

## plot ----
### auto ----
png(paste0("results/4_continuity/variance_decomposition_beta_SGDCs_noSeychelles_", comm, "_", comm_delin, "_bdmean.png"),
    width = 10, height = 5, units = 'in', res = 300)
par(oma = c(2,12,2,2),
    mar = c(1,1,1,1))
SGDC.decomp(SD = dist_merge[[metricSD]],
            GD = dist_merge[[metricGD]],
            FACTOR = dist_merge[,c("environment","seadist")])
par(oma = c(2,2,2,2),
    mar = c(1,1,1,1))
dev.off()



### perso plot ----
temp <- SGDC.decomp(SD = dist_merge[[metricSD]],
                    GD = dist_merge[[metricGD]],
                    FACTOR = dist_merge[,c("environment","seadist")])
temp <- 
  temp$Contribution %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("variable")

# rename variables
temp$variable <- gsub('Var\\(', '', temp$variable)
temp$variable <- gsub('Cov\\(', '', temp$variable)
temp$variable <- gsub('\\)', '', temp$variable)
temp$variable <- gsub('environment', 'Environment', temp$variable)
temp$variable <- gsub('seadist', 'Distance', temp$variable)
temp$variable <- gsub(',', ' ×', temp$variable)


temp$variable <- factor(temp$variable, levels = temp$variable[c(3,2,1,4,5)])

### {FIGURE 5} ####
ggplot(temp, aes(y = variable, x = Percent, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ffb734", "#9ECAE1", "#A1D99B", "grey", "grey20")) +
  theme_light() +
  ylab('') +
  theme(legend.position = "none",
        text = element_text(size = 15))
ggsave(paste0("results/4_continuity/_5_variance_decompositionPERSO_beta_SGDCs_noSeychelles_", comm, "_", comm_delin, "_bdmean.png"),
       width = 8, height = 4, dpi = 500)
ggsave(paste0("results/4_continuity/_5_variance_decompositionPERSO_beta_SGDCs_noSeychelles_", comm, "_", comm_delin, "_bdmean.pdf"),
       width = 8, height = 4)




# ---- 5. MRM ----

## MRM gd table ----
MRMgd <- tibble::tibble()
for (locations in c("allsites", "noSeychelles")){
  dist_merge <-
    read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))
  
  if(locations == "noSeychelles"){dist_merge <- mat.subset(dist_merge, "Seychelles")}
  
  for (metricGD in c("Fst", "GstPP.hed", "D.Jost")){
    
    statLM <- summary(lm(scale(dist_merge[[metricGD]]) ~
                           scale(dist_merge$seadist) +
                           scale(dist_merge$environment)))
    
    statMRM <- MRM(scale(dist_merge[[metricGD]]) ~
                     scale(dist_merge$seadist) +
                     scale(dist_merge$environment),
                   nperm = 9999)
    
    tempLM <-
      data.frame(statLM$coefficients[-1,1:2]) %>%
      rownames_to_column("explanatory_variable")
    
    tempMRM <-
      data.frame(pval = statMRM$coef[-1,2],
                 rsquared = statMRM$r.squared[1]) %>%
      rownames_to_column("explanatory_variable")
    # dplyr::rename(coefficients = "scale(dist_merge[[metricGD]])")
    # rbind(tibble(explanatory_variable = "R²",
    #              coefficients = statLM$adj.r.squared,
    #              pval = statMRM$r.squared[2])) %>%
    
    temp <-
      full_join(tempLM, tempMRM, by = "explanatory_variable") %>%
      cbind(locations = locations,
            response_variable = metricGD,
            .)
    
    MRMgd <-
      MRMgd %>%
      rbind(temp) %>%
      as_tibble()
  }
}

# rename variables 
MRMgd$explanatory_variable <- gsub('scale(dist_merge$', '', MRMgd$explanatory_variable, fixed = T)
MRMgd$explanatory_variable <- gsub(')', '', MRMgd$explanatory_variable, fixed = T)

# significance
MRMgd$signif <- ifelse(MRMgd$pval < 0.001, "***", 
                       ifelse(MRMgd$pval < 0.01, "**", 
                              ifelse(MRMgd$pval < 0.05, "*", "NS")))

# save
MRMgd %>%
  write_csv("results/4_continuity/MRM_beta_gd_envt_seadist_bdmean.csv")


## MRM sd table ----
MRMsd <- tibble::tibble()
for(comm_delin in comm_delin_list) {
  dist_merge <-
    read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))
  
  ### >>> Eupercaria ----
  colnames(dist_merge) <- gsub('Eupercaria/misc', 'Eupercaria', colnames(dist_merge))
  
  for (locations in c("allsites", "noSeychelles")){
    if(locations == "noSeychelles"){dist_merge <- mat.subset(dist_merge, "Seychelles")}
    
    for (comm in names_communities){
      
      for (metricSD in c("beta.jac", "beta.jtu", "beta.jne")){
        
        metricSDcomm <- paste(comm, metricSD, sep = ".")
        
        statLM <- summary(lm(scale(dist_merge[[metricSDcomm]]) ~
                               scale(dist_merge$seadist) +
                               scale(dist_merge$environment)))
        
        statMRM <- MRM(scale(dist_merge[[metricSDcomm]]) ~
                         scale(dist_merge$seadist) +
                         scale(dist_merge$environment),
                       nperm = 9999)
        
        tempLM <-
          data.frame(statLM$coefficients[-1,1:2]) %>%
          rownames_to_column("explanatory_variable")
        
        tempMRM <-
          data.frame(pval = statMRM$coef[-1,2],
                     rsquared = statMRM$r.squared[1]) %>%
          rownames_to_column("explanatory_variable")
        # dplyr::rename(coefficients = "scale(dist_merge[[metricSD]])")
        # rbind(tibble(explanatory_variable = "R²",
        #              coefficients = statLM$adj.r.squared,
        #              pval = statMRM$r.squared[2])) %>%
        
        temp <-
          full_join(tempLM, tempMRM, by = "explanatory_variable") %>%
          cbind(community_delineation = comm_delin,
                locations = locations,
                taxonomic_scale = comm,
                response_variable = metricSD,
                .)
        
        MRMsd <-
          MRMsd %>%
          rbind(temp) %>%
          as_tibble()
      }
    }
  }
}

# rename variables 
MRMsd$explanatory_variable <- gsub('scale(dist_merge$', '', MRMsd$explanatory_variable, fixed = T)
MRMsd$explanatory_variable <- gsub(')', '', MRMsd$explanatory_variable, fixed = T)

# significance
MRMsd$signif <- ifelse(MRMsd$pval < 0.001, "***", 
                       ifelse(MRMsd$pval < 0.01, "**", 
                              ifelse(MRMsd$pval < 0.05, "*", "NS")))
# save
MRMsd %>%
  write_csv("results/4_continuity/MRM_beta_sd_envt_seadist_bdmean.csv")



## DW plot ----
comm_delin <- comm_delin_list[1]

### load
MRMgd <- read_csv("results/4_continuity/MRM_beta_gd_envt_seadist_bdmean.csv")
MRMsd <- read_csv("results/4_continuity/MRM_beta_sd_envt_seadist_bdmean.csv")

MRMgd$response_variable <- factor(MRMgd$response_variable, levels = unique(MRMgd$response_variable))
MRMsd$response_variable <- factor(MRMsd$response_variable, levels = unique(MRMsd$response_variable))

MRMplot <-
  MRMgd %>% 
  bind_cols(taxonomic_scale = "Etelis coruscans") %>% 
  bind_cols(community_delineation = comm_delin) %>% 
  bind_rows(MRMsd) %>% 
  dplyr::filter(!response_variable %in% c("GstPP.hed", "D.Jost", "beta.jne")) %>% 
  dplyr::filter(community_delineation == comm_delin)

MRMplot$taxonomic_scale <- 
  factor(MRMplot$taxonomic_scale, levels = c("Etelis coruscans", names_communities))

MRMplot <- 
  MRMplot %>% 
  mutate(locations = replace(locations, locations == "allsites", "all sites")) %>% 
  mutate(locations = replace(locations, locations == "noSeychelles", "without Seychelles"))

color_response <- c(Fst =  "#F69C73FF",
                    # D.Jost = "#03051AFF",
                    beta.jac = "#60CEACFF",
                    beta.jtu = "#395D9CFF")


### {FIGURE S7} ####
ggplot(MRMplot, aes(x=Estimate, y=explanatory_variable, 
                     group=response_variable, color=response_variable, shape=response_variable)) +
  geom_vline(xintercept = 0,linetype=2) +
  geom_pointrange(aes(xmin=Estimate - 1.96*Std..Error,
                      xmax=Estimate + 1.96*Std..Error),
                  position=position_dodge(width = -0.4)) +
  scale_color_manual('', values = color_response)+
  facet_grid(rows = vars(taxonomic_scale), 
             cols = vars(locations)) +
  xlim(c(-1,1)) + 
  ylab("") +
  labs(shape = "") +
  theme_light() +
  theme(strip.background = element_rect(color = "grey"),
        legend.position = "top")
  # theme(axis.text.y = element_text(color = color_explanatory))

ggsave(paste0("results/4_continuity/_S7_DWplotPERSO_MRM_beta_combined_envt_seadist.png"),
         height = 7, width = 6, units = 'in', dpi = 500)





# 6. SGDCs residuals ----
# multi.mantel: same values as MRM ! but allows to take residuals

comm_delin = comm_delin_list[1]

dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))
dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))

dist_mat <- lapply(dist_mat, function(x) mat.subset(x, "Seychelles"))
dist_merge <- mat.subset(dist_merge, "Seychelles")

metricGD = "Fst"
metricSD = "beta.jtu"
metricXLAB <- gsub('beta.', 'β', metricSD)

gglist <- list()
for (comm in names_communities){
  
  # comm = names_communities[2]
  metricSDcomm = paste(comm, metricSD, sep = ".")
  
  
  ## MRM ----
  modelGD <-
    phytools::multi.mantel((dist_mat[[metricGD]]),
                           list(dist_mat$environment,
                                dist_mat$seadist),
                           nperm = 1000)
  
  modelSD <-
    phytools::multi.mantel((dist_mat[[metricSDcomm]]),
                           list(dist_mat$environment,
                                dist_mat$seadist),
                           nperm = 1000)
  
  stats_res <- 
    mantel(modelGD$residuals, modelSD$residuals, permutations = 9999)
  
  
  
  ## LM ----
  modelGD <-
    lm((dist_merge[[metricGD]]) ~ 
         dist_merge$environment +
         dist_merge$seadist)
  
  modelSD <-
    lm((dist_merge[[metricSDcomm]]) ~
         dist_merge$environment +
         dist_merge$seadist)
  
  summary(lm(modelGD$residuals ~ modelSD$residuals))
  
  
  df <- data.frame(SDres = modelSD$residuals, GDres = modelGD$residuals)
  
  if(stats_res$signif <= 0.05){
    gglist[[comm]] <-
      ggplot(df, aes(x = SDres, y = GDres)) +
      geom_smooth(method = "lm", color = "grey50", fill = "grey80") +
      geom_point(shape = 19, alpha = 0.8) +
      xlab(paste0("residuals(", metricXLAB, ") (", comm, ")")) +
      ylab(paste0("residuals(", metricGD, ")")) +
      annotate('text',
               x = min(df$SDres), y = 1.5*max(df$GDres),
               hjust = 0, vjust = 1,
               label = paste0("r Mantel = ", round(stats_res$statistic, 4), "\np = ", stats_res$signif)) +
      theme_light()
  }
  
  if(stats_res$signif > 0.05){
    gglist[[comm]] <-
      ggplot(df, aes(x = SDres, y = GDres)) +
      geom_point(shape = 19, alpha = 0.8) +
      xlab(paste0("residuals(", metricXLAB, ") (", comm, ")")) +
      ylab(paste0("residuals(", metricGD, ")")) +
      annotate('text',
               x = min(df$SDres), y = 1.5*max(df$GDres),
               hjust = 0, vjust = 1,
               label = paste0("r Mantel = ", round(stats_res$statistic, 4), "\np = ", stats_res$signif)) +
      theme_light()
  }
}

### {FIGURE S10} ####
gg_grob <-
  patchwork::wrap_plots(gglist, tag_level = "new") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
plot(gg_grob)

ggsave(width = 8, height = 7, dpi = 500,
       filename = paste0("results/4_continuity/_S10_beta_SGDCs_plot_significant_RESIDUALS.png"))
ggsave(width = 8, height = 7,
       filename = paste0("results/4_continuity/_S10_beta_SGDCs_plot_significant_RESIDUALS.pdf"))









# ---- 7. dbRDA (geo + env) ----
## dbRDA GD table ----

varpartGD <- tibble::tibble()
dbRDAgd <- tibble::tibble()

for (locations in c("All sites", "Without Seychelles")){
  dist_merge <-
    read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"), col_types = cols())
  
  dist_mat <-
    readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))
  
  envt_site <- 
    read_csv("intermediate/3_distance_metrics/bio_oracle_variables_bdmean_stationbuffer_siteaverage.csv", col_types = cols()) 
  
  coord_sites <-
    data_sites %>%
    dplyr::select(-number_samples) %>% 
    shift.lon()
  
  if(locations == "Without Seychelles"){
    dist_merge <- mat.subset(dist_merge, "Seychelles")
    dist_mat <- lapply(dist_mat, function(x) mat.subset(x, "Seychelles"))
    envt_site <- mat.subset(envt_site, "Seychelles")
    coord_sites <- mat.subset(coord_sites, "Seychelles") %>% column_to_rownames("site")
  }
  
  ### a. dbMEM on seadist
  dbMEM_seadist <- adespatial::dbmem(dist_mat$seadist, store.listw = TRUE) # MEM.autocor = "non-null",
  # dbMEM_seadist <- adespatial::dbmem(coord_sites, store.listw = TRUE) # MEM.autocor = "non-null",
  
  ## keep fist 2 axis, explains ~95% of variation
  # PCA_seadist <- prcomp(dist_mat$seadist)
  # factoextra::fviz_eig(PCA_seadist)
  # sum((PCA_seadist$sdev^2/sum(PCA_seadist$sdev^2))[1:3])
  dbMEM_seadist <- dbMEM_seadist[,1:2] 
  

  ### b. environment PCA (X)
  ## scale
  envt_site_scaled <-
    envt_site %>% 
    column_to_rownames("site") %>% 
    scale(center = T, scale = T) %>% 
    as.data.frame()
  
  ## PCA
  PCA_envt <- prcomp(envt_site_scaled)
  PCA_envtaxis <- as.data.frame(PCA_envt$x) # The principal components can be found in the $x matrix

  ## keep fist 2 axis, explains ~95% of variation
  # sum((PCA_envt$sdev^2/sum(PCA_envt$sdev^2))[1:2])
  # factoextra::fviz_eig(PCA_envt)
  PCA_envtaxis <- PCA_envtaxis[,1:2]
  
  
  for (metricGD in c("Fst", "GstPP.hed", "D.Jost")){ #

    ### a. dbMEM on seadist select variables
    dbrda_GD0 <- vegan::capscale(dist_mat[[metricGD]] ~ 1, data = dbMEM_seadist) # null model, only intercept
    dbrda_GDgeo <- vegan::capscale(dist_mat[[metricGD]] ~ ., data = dbMEM_seadist)
    selGEO <- ordistep(dbrda_GD0, scope = formula(dbrda_GDgeo), 
                       direction="both", trace = FALSE)
    
    MEMstokeep <- 
      rownames(selGEO$anova) %>% 
      grep(pattern = "\\+ ", value = TRUE) %>% 
      gsub(pattern = "\\+ ", replacement = "")
    
    dbMEM_seadistsel <-
      dbMEM_seadist %>% 
      as.data.frame() %>% 
      dplyr::select(all_of(MEMstokeep))
    
    
    ### best model table
    X = cbind(dbMEM_seadistsel, PCA_envtaxis)
    dbrda_GD0 <- vegan::capscale(dist_mat[[metricGD]] ~ 1, data = X) # null model, only intercept
    dbrda_GDall <- vegan::capscale(dist_mat[[metricGD]] ~ ., data = X)
    
    dbrda_GDsel <- ordistep(dbrda_GD0, scope = formula(dbrda_GDall), 
                            direction="both", trace = FALSE) 
    # dbrda_GDsel$anova
    
    if(!is.null(dbrda_GDsel$anova)){
      dbRDAgd <-
        dbRDAgd %>%
        rbind(tibble(location = locations,
                     Model = as.character(dbrda_GDsel$call)[2],
                     R2 = RsquareAdj(dbrda_GDsel)$r.squared, 
                     adjR2 = RsquareAdj(dbrda_GDsel)$adj.r.squared,
                     F = anova.cca(dbrda_GDsel)[1,3],
                     pval = anova.cca(dbrda_GDsel)[1,4])) %>% 
        mutate(Model = gsub("dist_mat\\[\\[metricGD\\]\\]", metricGD, Model))
    }
    
    if(is.null(dbrda_GDsel$anova)){
      dbRDAgd <-
        dbRDAgd %>%
        rbind(tibble(location = locations,
                     Model = metricGD,
                     R2 = NA, 
                     adjR2 = NA,
                     F = NA,
                     pval = NA))
    }
    
    
    ### varpart
    # vpGD <- varpart(dist_mat[[metricGD]], dbMEM_seadist, PCA_envtaxis)
    # vpGD
    # showvarparts(2, bg = c("#9ECAE1", "#A1D99B"))
    # plot(vpGD, bg = c("#9ECAE1", "#A1D99B"), Xnames = c("Distance", "Environment"),
    #      main = metricGD)
    

    ## Varpart computatio by hand : sligthly different from varpart plot...
    RDAabc <- vegan::capscale(dist_mat[[metricGD]] ~ ., data = cbind(dbMEM_seadist, PCA_envtaxis))
    RDAa <- vegan::capscale(dist_mat[[metricGD]] ~ MEM1 + MEM2 + Condition(PC1 + PC2), data = cbind(dbMEM_seadist, PCA_envtaxis))
    RDAb <- vegan::capscale(dist_mat[[metricGD]] ~ PC1 + PC2 + Condition(MEM1 + MEM2), data = cbind(dbMEM_seadist, PCA_envtaxis))
    
    a <- RsquareAdj(RDAa)$adj.r.squared
    b <- RsquareAdj(RDAb)$adj.r.squared
    c <- RsquareAdj(RDAabc)$adj.r.squared - a - b
    
    varpartGD <- 
      varpartGD %>% 
      rbind(tibble(location = locations,
                   response_variable = metricGD,
                   explanatory_variable = c("Distance", "Environment", "Shared", "Residuals"),
                   variation_explained = c(a, b, c, 1-a-b-c),
                   pval = c(anova(RDAa)[1,4], anova(RDAb)[1,4], NA, NA))
      )
  }
}


# significance
dbRDAgd$signif <- ifelse(dbRDAgd$pval < 0.001, "***", 
                         ifelse(dbRDAgd$pval < 0.01, "**", 
                                ifelse(dbRDAgd$pval < 0.05, "*", "NS")))

varpartGD$signif <- ifelse(varpartGD$pval < 0.001, "***", 
                         ifelse(varpartGD$pval < 0.01, "**", 
                                ifelse(varpartGD$pval < 0.05, "*", "NS")))
dbRDAgd # different runs of ordistep can give different results
varpartGD

varpartGDwide <- 
  varpartGD %>% 
  mutate(variation_explained = paste(round(variation_explained, 4), signif)) %>% 
  mutate(variation_explained = gsub(" NA", "", variation_explained)) %>% 
  dplyr::select(-c(pval, signif)) %>% 
  pivot_wider(names_from = explanatory_variable, values_from = variation_explained)

# save
dbRDAgd %>% write_csv("results/4_continuity/dbRDA_models_beta_gd.csv")
varpartGDwide %>% write_csv("results/4_continuity/dbRDA_varpart_beta_gd.csv")





## dbRDA SD table ----

varpartSD <- tibble::tibble()
dbRDAsd <- tibble::tibble()
for (locations in c("All sites", "Without Seychelles")){
  dist_merge <-
    read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"), col_types = cols())
  
  dist_mat <-
    readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))
  
  envt_site <- 
    read_csv("intermediate/3_distance_metrics/bio_oracle_variables_bdmean_stationbuffer_siteaverage.csv", col_types = cols()) 
  
  if(locations == "Without Seychelles"){
    dist_merge <- mat.subset(dist_merge, "Seychelles")
    dist_mat <- lapply(dist_mat, function(x) mat.subset(x, "Seychelles"))
    envt_site <- mat.subset(envt_site, "Seychelles")
  }
  
  ### a. dbMEM on seadist
  dbMEM_seadist <- adespatial::dbmem(dist_mat$seadist, store.listw = TRUE) # MEM.autocor = "non-null",
  # dbMEM_seadist <- adespatial::dbmem(coord_sites, store.listw = TRUE) # MEM.autocor = "non-null",
  
  ## keep fist 2 axis, explains ~95% of variation
  # PCA_seadist <- prcomp(dist_mat$seadist)
  # factoextra::fviz_eig(PCA_seadist)
  # sum((PCA_seadist$sdev^2/sum(PCA_seadist$sdev^2))[1:2])
  dbMEM_seadist <- dbMEM_seadist[, 1:2] 
  
  
  ### b. environment PCA (X)
  ## scale
  envt_site_scaled <-
    envt_site %>% 
    column_to_rownames("site") %>% 
    scale(center = T, scale = T) %>% 
    as.data.frame()
  
  ## PCA
  PCA_envt <- prcomp(envt_site_scaled)
  PCA_envtaxis <- as.data.frame(PCA_envt$x) # The principal components can be found in the $x matrix
  # sum((PCA_envt$sdev^2/sum(PCA_envt$sdev^2))[1:2])
  # factoextra::fviz_eig(PCA_envt)
  
  ## keep fist 2 axis, explains ~95% of variation
  PCA_envtaxis <- PCA_envtaxis[, 1:2]
  
  
  for (comm in names_communities){
    
    for (metricSD in c("beta.jac", "beta.jtu")){ #"beta.jne"
      
      metricSDcomm <- paste(comm, metricSD, sep = ".")
      
      ### a. dbMEM on seadist select variables
      dbrda_SD0 <- vegan::capscale(dist_mat[[metricSDcomm]] ~ 1, data = dbMEM_seadist) # null model, only intercept
      dbrda_SDgeo <- vegan::capscale(dist_mat[[metricSDcomm]] ~ ., data = dbMEM_seadist)
      selGEO <- ordistep(dbrda_SD0, scope = formula(dbrda_SDgeo), 
                         direction="both", trace = FALSE)
      
      MEMstokeep <- 
        rownames(selGEO$anova) %>% 
        grep(pattern = "\\+ ", value = TRUE) %>% 
        gsub(pattern = "\\+ ", replacement = "")
      
      dbMEM_seadistsel <-
        dbMEM_seadist %>% 
        as.data.frame() %>% 
        dplyr::select(all_of(MEMstokeep))
      
      
      
      ### best model table
      X = cbind(dbMEM_seadistsel, PCA_envtaxis)
      dbrda_SD0 <- vegan::capscale(dist_mat[[metricSDcomm]] ~ 1, data = X) # null model, only intercept
      dbrda_SDall <- vegan::capscale(dist_mat[[metricSDcomm]] ~ ., data = X)
      
      dbrda_SDsel <- ordistep(dbrda_SD0, scope = formula(dbrda_SDall), 
                              direction="both", trace = FALSE) 
      # dbrda_SDsel$anova
      
      if(!is.null(dbrda_SDsel$anova)){
        dbRDAsd <-
          dbRDAsd %>% 
          rbind(tibble(location = locations,
                       Model = as.character(dbrda_SDsel$call)[2],
                       R2 = RsquareAdj(dbrda_SDsel)$r.squared, 
                       adjR2 = RsquareAdj(dbrda_SDsel)$adj.r.squared,
                       F = anova.cca(dbrda_SDsel)[1,3],
                       pval = anova.cca(dbrda_SDsel)[1,4])) %>% 
          mutate(Model = gsub("dist_mat\\[\\[metricSDcomm\\]\\]", metricSDcomm, Model))
      }
      
      if(is.null(dbrda_SDsel$anova)){
        dbRDAsd <-
          dbRDAsd %>% 
          rbind(tibble(location = locations,
                       Model = metricSDcomm,
                       R2 = NA, 
                       adjR2 = NA,
                       F = NA,
                       pval = NA))
      }
        
        ### varpart
        # vpSD <- varpart(dist_mat[[metricSDcomm]], dbMEM_seadist, PCA_envtaxis)
        # vpSD
        # showvarparts(2, bg = c("#9ECAE1", "#A1D99B"))
        # plot(vpSD, bg = c("#9ECAE1", "#A1D99B"), Xnames = c("Distance", "Environment"),
        #      main = metricSDcomm)
        
        
        ## Varpart computatio by hand : sligthly different from varpart plot...
        RDAabc <- vegan::capscale(dist_mat[[metricSDcomm]] ~ ., data = cbind(dbMEM_seadist, PCA_envtaxis))
        RDAa <- vegan::capscale(dist_mat[[metricSDcomm]] ~ MEM1 + MEM2 + Condition(PC1 + PC2), data = cbind(dbMEM_seadist, PCA_envtaxis))
        RDAb <- vegan::capscale(dist_mat[[metricSDcomm]] ~ PC1 + PC2 + Condition(MEM1 + MEM2), data = cbind(dbMEM_seadist, PCA_envtaxis))
        
        a <- RsquareAdj(RDAa)$adj.r.squared
        b <- RsquareAdj(RDAb)$adj.r.squared
        c <- RsquareAdj(RDAabc)$adj.r.squared - a - b
        
        varpartSD <- 
          varpartSD %>% 
          rbind(tibble(location = locations,
                       response_variable = metricSDcomm,
                       explanatory_variable = c("Distance", "Environment", "Shared", "Residuals"),
                       variation_explained = c(a, b, c, 1-a-b-c),
                       pval = c(anova(RDAa)[1,4], anova(RDAb)[1,4], NA, NA))
          )
        
    }
  }
}

# significance
dbRDAsd$signif <- ifelse(dbRDAsd$pval < 0.001, "***", 
                         ifelse(dbRDAsd$pval < 0.01, "**", 
                                ifelse(dbRDAsd$pval < 0.05, "*", "NS")))

varpartSD$signif <- ifelse(varpartSD$pval < 0.001, "***", 
                           ifelse(varpartSD$pval < 0.01, "**", 
                                  ifelse(varpartSD$pval < 0.05, "*", "NS")))
dbRDAsd # different runs of ordistep can give different results

varpartSDwide <- 
  varpartSD %>% 
  mutate(variation_explained = paste(round(variation_explained, 4), signif)) %>% 
  mutate(variation_explained = gsub(" NA", "", variation_explained)) %>% 
  dplyr::select(-c(pval, signif)) %>% 
  pivot_wider(names_from = explanatory_variable, values_from = variation_explained)
varpartSDwide

# save
dbRDAsd %>% write_csv("results/4_continuity/dbRDA_models_beta_sd.csv")
varpartSDwide %>% write_csv("results/4_continuity/dbRDA_varpart_beta_sd.csv")










# 
# ---- *** DRAFTS ----
## *** SGDCs by distance phylo ----
# 
# SGDCs <- data.frame(community = names_communities)
# 
# 
# temp <- lapply(SGDCs$community, FUN = function(x){stat_Mantel[[x]]$statistic})
# names(temp) <- names_communities
# 
# SGDCs <- 
#   temp %>% 
#   as.data.frame() %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var = "community") %>% 
#   dplyr::rename(rMantel = V1)
# 
# SGDCs$phylodist <- 
#   gsub(pattern = "phylodist_", replacement = "", SGDCs$community)
# 
# plot(SGDCs$phylodist, SGDCs$rMantel)
# summary(lm(SGDCs$phylodist ~ SGDCs$rMantel))
# 
# 
# 
# 
## *** DW plot automatic ----
# 
# ### GD
# models_gd <- list()
# for (metricGD in c("Fst", "GstPP.hed", "D.Jost")){
#   models_gd[[metricGD]] <- lm(scale(dist_merge[[metricGD]]) ~
#                                 scale(dist_merge$environment) +
#                                 scale(dist_merge$seadist))
#   
# }
# dwplot(models_gd,
#        vline = geom_vline(xintercept = 0,
#                           colour = "grey60",
#                           linetype = 2))
# 
# ggsave("results/4_continuity/DWplot_MRM_beta_gd_envt_seadist_noSeychelles.png",
#        height = 5, width = 8, units = 'in', dpi = 300)
# 
# ### SD
# models_sd <- list()
# for (metric in c("beta.jac", "beta.jtu", "beta.jne")){
#   for (comm in names_communities){
#     
#     metricSD <- paste(comm, metric, sep = ".")
#     
#     models_sd[[metricSD]] <- lm(scale(dist_merge[[metricSD]]) ~
#                                   scale(dist_merge$environment) +
#                                   scale(dist_merge$seadist))
#   }
# }
# dwplot(models_sd,
#        vline = geom_vline(xintercept = 0,
#                           colour = "grey60",
#                           linetype = 2))
# ggsave(paste0("results/4_continuity/DWplot_MRM_beta_sd_envt_seadist_noSeychelles_", comm_delin, ".png"),
#        height = 5, width = 8, units = 'in', dpi = 300)
#
## *** cor explanatory variables ----
# model <- cor(dist_merge[c("geodist",
#                "seadist",
#                "leastcost",
#                "environment")])
# corrplot::corrplot(model)
# 
# cor(dist_merge[["seadist"]], dist_merge[["environment"]])
# 
# 
# 
# 
# 
## ---- *** old plot IBD ----
# 
# # # MRM
# # sMRM_IBDsd <- MRM(dist_mat[[metricSD]] ~ dist_mat[[metricDIST]], nperm = 9999)
# # sMRM_IBDgd <- MRM(dist_mat[[metricGD]] ~ dist_mat[[metricDIST]], nperm = 9999)
# # sMRM_SGDC <- MRM(dist_mat[[metricSD]] ~ dist_mat[[metricGD]], nperm = 9999)
# # # MRM(dist_merge[[metricSD]] ~ dist_merge[[metricDIST]], nperm = 9999) # same with distance matrix or long df
# 
# # Mantel
# sMantel_IBDsd <- vegan::mantel(dist_mat[[metricSD]], dist_mat[[metricDIST]], permutations = 9999)
# sMantel_IBDgd <- vegan::mantel(dist_mat[[metricGD]], dist_mat[[metricDIST]], permutations = 9999)
# sMantel_SGDC <- vegan::mantel(dist_mat[[metricSD]], dist_mat[[metricGD]], permutations = 9999)
# 
# # # distance decay models
# # sDecay_IBDsd <- decay.model(dist_mat[[metricSD]], dist_mat[[metricDIST]], y.type="dissim", model.type="pow", perm=999)
# # sDecay_IBDgd <- decay.model(dist_mat[[metricGD]], dist_mat[[metricDIST]], y.type="dissim", model.type="pow", perm=999)
# #
# # plot(dist_mat[[metricDIST]], dist_mat[[metricSD]], pch = 16)
# # plot.decay(sDecay_IBDsd, col="green", remove.dots=T, add=T, lwd=2)
# #
# # plot(dist_mat[[metricDIST]], dist_mat[[metricGD]], pch = 16)
# # plot.decay(sDecay_IBDgd, col="green", remove.dots=T, add=T, lwd=2)
# 
# 
# 
# # add Seychelles color
# dist_merge$Seychelles <- "No"
# dist_merge[grep("Seychelles", dist_merge$site),]$Seychelles <- "Yes"
# 
# # species IBD
# ggSD <-
#   ggplot(dist_merge) +
#   geom_point(aes(.data[[metricDIST]], .data[[metricSD]], color = Seychelles), show.legend = F) +
#   scale_color_brewer(palette="Dark2") +
#   annotate('text',
#            x=min(dist_merge[[metricDIST]]), y=max(dist_merge[[metricSD]]),
#            hjust = 0, vjust = 1,
#            label=paste0("r Mantel = ", sprintf("%.3f", sMantel_IBDsd$statistic), ", p = ", sMantel_IBDsd$signif)) +
#   # "\nR² MRM = ", round(sMRM_IBDsd$r.squared[["R2"]], 4), ", p = ", round(sMRM_IBDsd$r.squared[["pval"]], 5))) +
#   theme_light()
# 
# # genetic IBD
# ggGD <-
#   ggplot(dist_merge) +
#   geom_point(aes(.data[[metricDIST]], .data[[metricGD]], color = Seychelles), show.legend = F) +
#   scale_color_brewer(palette="Dark2") +
#   annotate('text',
#            x=min(dist_merge[[metricDIST]]), y=max(dist_merge[[metricGD]]),
#            hjust = 0, vjust = 1,
#            label=paste0("r Mantel = ", sprintf("%.3f", sMantel_IBDgd$statistic), ", p = ", sMantel_IBDgd$signif)) +
#   # "\nr MRM = ", round(sMRM_IBDgd$r.squared[["R2"]], 4), ", p = ", round(sMRM_IBDgd$r.squared[["pval"]], 5))) +
#   theme_light()
# 
# 
# # SGDCs
# ggSGDCs <-
#   ggplot(dist_merge) +
#   geom_point(aes(.data[[metricSD]], .data[[metricGD]], color = Seychelles)) +
#   scale_color_brewer(palette="Dark2") +
#   annotate('text',
#            x=min(dist_merge[[metricSD]]), y=max(dist_merge[[metricGD]]),
#            hjust = 0, vjust = 1,
#            label=paste0("r Mantel = ", sprintf("%.3f", sMantel_SGDC$statistic), ", p = ", sMantel_SGDC$signif)) +
#   # "\nr MRM = ", round(sMRM_SGDC$r.squared[["R2"]], 4), ", p = ", round(sMRM_SGDC$r.squared[["pval"]], 5))) +
#   theme_light()
# 
# 
# # merge plots
# ggSD + ggGD + ggSGDCs + #+ plot_annotation(title = comm)
#   plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
# ggsave(width = 14, height = 4,
#        filename = paste0("results/4_continuity/IBD_beta_SGDC_noSeychelles_", metricSD, "_", metricGD, "_", metricDIST, "_", comm_delin, ".png"))
# 


## *** GD ----
# ggplot(MRMgd, aes(x=Estimate, y=explanatory_variable, group=response_variable, color=response_variable)) +
#   geom_vline(xintercept = 0,linetype=2) +
#   geom_pointrange(aes(xmin=Estimate - 1.96*Std..Error,
#                       xmax=Estimate + 1.96*Std..Error),
#                   position=position_dodge(width = -0.4)) +
#   scale_color_viridis(option="rocket", discrete=T, end = 0.8)+
#   facet_wrap(~locations) +
#   xlim(c(-1,1)) +
#   theme_light()
# 
# ggsave(paste0("results/4_continuity/DWplotPERSO_MRM_beta_gd_envt_seadist_bdmean.png"),
#        height = 5, width = 12, units = 'in', dpi = 300)
# 
# 
## *** SD ----
# # all community delineations
# ggplot(MRMsd, aes(x=Estimate,
#                   y=explanatory_variable,
#                   shape=community_delineation,
#                   group=interaction(community_delineation, response_variable),
#                   color=response_variable)) +
#   geom_vline(xintercept = 0, linetype=2) +
#   geom_pointrange(aes(xmin=Estimate - 1.96*Std..Error,
#                     xmax=Estimate + 1.96*Std..Error),
#                     position=position_dodge(width=-1)) +
#   scale_color_viridis(option="mako", discrete=T, end = 0.8) +
#   facet_grid(vars(taxonomic_scale), vars(locations), scale="free") +
#   xlim(c(-1,1)) +
#   theme_light()
# 
# ggsave(paste0("results/4_continuity/DWplotPERSO_MRM_beta_sd_envt_seadist_bdmean.png"),
#        height = 10, width = 12, units = 'in', dpi = 300)
# 
# 
# # metric by metric
# for(metricSD in unique(MRMsd$response_variable)){
#   MRMsub <-
#     MRMsd %>%
#     filter(response_variable == metricSD)
# 
#   ggplot(MRMsub, aes(x=Estimate,
#                     y=explanatory_variable,
#                     shape=community_delineation,
#                     color=community_delineation)) +
#     geom_vline(xintercept = 0, linetype=2) +
#     geom_pointrange(aes(xmin=Estimate - 1.96*Std..Error,
#                         xmax=Estimate + 1.96*Std..Error),
#                     position=position_dodge(width=-0.5)) +
#     scale_color_viridis(option="cividis", discrete=T, end = 0.8) +
#     facet_grid(vars(taxonomic_scale), vars(locations), scale="free") +
#     xlim(c(-1,1)) +
#     theme_light() +
#     ggtitle(metricSD)
# 
#   ggsave(paste0("results/4_continuity/DWplotPERSO_MRM_beta_sd_envt_seadist_", metricSD, "_bdmean.png"),
#          height = 8, width = 12, units = 'in', dpi = 300)
# 
# }
# 
# 
# 
# # only taxonomy
# comm_delin <- "taxonomy"
# 
# MRMsub <-
#   MRMsd %>%
#   filter(community_delineation == comm_delin)
# 
# MRMsub$response_variable <- 
#   factor(MRMsub$response_variable,
#          level = c("beta.jne", "beta.jac", "beta.jtu"))
# 
# ggplot(MRMsub, aes(x=Estimate,
#                   y=explanatory_variable,
#                   group=interaction(response_variable),
#                   color=response_variable)) +
#   geom_vline(xintercept = 0, linetype=2) +
#   geom_pointrange(aes(xmin=Estimate - 1.96*Std..Error,
#                       xmax=Estimate + 1.96*Std..Error),
#                   position=position_dodge(width=-0.4)) +
#   scale_color_manual(values = c("#0B0405FF", "#60CEACFF", "#395D9CFF"))+
#   # scale_color_viridis(option="mako", discrete=T, end = 0.8) +
#   facet_grid(vars(taxonomic_scale), vars(locations), scale="free") +
#   xlim(c(-1,1)) +
#   theme_light()
# 
# ggsave(paste0("results/4_continuity/DWplotPERSO2_MRM_beta_sd_envt_seadist_bdmean_", comm_delin, ".png"),
#        height = 8, width = 10, units = 'in', dpi = 300)
# 
# 
# 
# 
# 
# 
# 
## *** Varapart computation ----
## METHOD 1 computation
# RDAabc <- vegan::capscale(dist_mat[[metricGD]] ~ ., data = cbind(dbMEM_seadist, PCA_envtaxis))
# RDAac <- vegan::capscale(dist_mat[[metricGD]] ~ ., data = dbMEM_seadist)
# RDAbc <- vegan::capscale(dist_mat[[metricGD]] ~ ., data = PCA_envtaxis)

# abc <- RsquareAdj(RDAabc)$adj.r.squared # = a+b+c 
# ac <- RsquareAdj(RDAac)$adj.r.squared   # = a+c seadist
# bc <- RsquareAdj(RDAbc)$adj.r.squared   # = b+c environment
# c <- ac + bc - abc  # covariance

# a <- ac - c         # seadist
# b <- bc - c         # environment
