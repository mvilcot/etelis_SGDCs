# ---- Setup ----
## parameters ----

names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")

sites = "noSeychelles"

comm_delin <- "taxonomy"


## load ----
dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))

envt_site <- 
  read_csv("intermediate/3_distance_metrics/bio_oracle_variables_bdmean_stationbuffer_siteaverage.csv") # %>%




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



# ---- 1. IBD and IBE ----
## a. IBD ----
metricGD <- "Fst"
metricSDs <- paste(names_communities, "beta.jtu", sep = ".")
metricDIST <- 'seadist' #'environment'

gg_list <- list()
IBD_table <- tibble()
i=1

library(lmPerm)

for (metricDIV in c(metricGD, metricSDs)){
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
    
    ## Procruste
    sProtest <- vegan::protest(dist_matSUB[[metricDIV]], dist_matSUB[[metricDIST]], permutations = 9999)
    sProtest
    
    ## PLOT
    # color
    color_x <- "#3182BD"
    if(metricDIV %in% metricGD) color_y <- "#F69C73FF" # rocket
    if(metricDIV %in% metricSDs) color_y <- "#395D9CFF" # mako
    
    # proper beta and families
    metricLAB <- gsub('beta.', 'β', metricDIV)
    if(metricDIV %in% metricSDs){
      temp <- str_split_1(metricLAB, '[.]')
      metricLAB <- paste0(temp[2], ' (', temp[1], ')')
    }
    
    # plot linear model if significant
    if(sMantel$signif < 0.05){
      gg_list[[i]] <-
        ggplot(dist_mergeSUB, aes(.data[[metricDIST]], .data[[metricDIV]])) + #, color = Seychelles
        geom_smooth(method = "lm", color = "grey50", se = F) +
        geom_point(aes(shape = Seychelles), show.legend = F, color = color_y) +
        scale_shape_manual(values = c(19, 1)) +
        annotate('text',
                 x=min(dist_mergeSUB[[metricDIST]]), y=max(dist_mergeSUB[[metricDIV]]),
                 hjust = 0, vjust = 1,
                 label=paste0("r = ", sprintf("%.3f", sMantel$statistic), "\np = ", sMantel$signif)) +
        theme_light() +
        ylab(metricLAB) + xlab('Distance by sea') +
        theme(
          axis.title.x = element_text(colour = color_x),
          # axis.title.y = element_text(colour = color_y)
        )
    }
    
    # no linear model if not significant
    if(sMantel$signif > 0.05 | (locations == "allsites" & metricDIV == "Fst")){
      gg_list[[i]] <-
        ggplot(dist_mergeSUB, aes(.data[[metricDIST]], .data[[metricDIV]], shape = Seychelles)) + #
        geom_point(aes(shape = Seychelles), show.legend = F, color = color_y) +
        scale_shape_manual(values = c(19, 1)) +
        annotate('text',
                 x=min(dist_mergeSUB[[metricDIST]]), y=max(dist_mergeSUB[[metricDIV]]),
                 hjust = 0, vjust = 1,
                 label=paste0("r = ", sprintf("%.3f", sMantel$statistic), "\np = ", sMantel$signif)) +
        theme_light() +
        ylab(metricLAB) + xlab('Distance by sea') +
        theme(
          axis.title.x = element_text(colour = color_x),
          # axis.title.y = element_text(colour = color_y)
        )
    }
    
    IBD_table <- 
      IBD_table %>% 
      rbind(tibble(locations = locations,
                   metric = metricDIV,
                   isolation_by = metricDIST,
                   Mantel_r = sMantel$statistic,
                   Mantel_pval = sMantel$signif,
                   Procruste_r = sqrt(1-sProtest$ss),
                   Procruste_pval = sProtest$signif))
    
    i = i+1
  }
}


### {FIGURE S9} ####
gg_grob <-
  patchwork::wrap_plots(gg_list, ncol = 2, tag_level = "new") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
plot(gg_grob)

ggsave(gg_grob, width = 10, height = 15, dpi = 500,
       filename = paste0("results/4_continuity/_S7_isolation_by_", metricDIST, "_", comm_delin, "_2.png"))


### {TABLE S5} ####
IBD_table$Mantel_signif <- ifelse(IBD_table$Mantel_pval < 0.001, "***", 
                                  ifelse(IBD_table$Mantel_pval < 0.01, "**", 
                                         ifelse(IBD_table$Mantel_pval < 0.05, "*", "NS")))
IBD_table$Procruste_signif <- ifelse(IBD_table$Procruste_pval < 0.001, "***", 
                                     ifelse(IBD_table$Procruste_pval < 0.01, "**", 
                                            ifelse(IBD_table$Procruste_pval < 0.05, "*", "NS")))
IBD_table <-
  IBD_table %>%
  mutate(Mantel_r = paste(round(Mantel_r, 3), Mantel_signif)) %>%
  mutate(Procruste_r = paste(round(Procruste_r, 3), Procruste_signif)) %>%
  mutate(Mantel_r = gsub(" NS", "", Mantel_r)) %>%
  mutate(Procruste_r = gsub(" NS", "", Procruste_r)) %>%
  dplyr::select(-c(Mantel_pval, Procruste_pval, Mantel_signif, Procruste_signif)) %>%
  dplyr::arrange(locations)

write_csv(IBD_table, paste0("results/4_continuity/_Sx_isolation_by_", metricDIST, "_", comm_delin, ".csv"))



## b. IBE ----
metricGD <- "Fst"
metricSDs <- paste(names_communities, "beta.jtu", sep = ".")
metricDIST <- 'environment'

gg_list <- list()
IBE_table <- tibble()
i=1

for (metricDIV in c(metricGD, metricSDs)){
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
    
    ## Procruste
    sProtest <- vegan::protest(dist_matSUB[[metricDIV]], dist_matSUB[[metricDIST]], permutations = 9999)
    sProtest
    
    ## PLOT
    # color
    color_x <- "#31A354"
    if(metricDIV %in% metricGD) color_y <- "#F69C73FF" # rocket
    if(metricDIV %in% metricSDs) color_y <- "#395D9CFF" # mako
    
    # proper beta and families
    metricLAB <- gsub('beta.', 'β', metricDIV)
    if(metricDIV %in% metricSDs){
      temp <- str_split_1(metricLAB, '[.]')
      metricLAB <- paste0(temp[2], ' (', temp[1], ')')
    }
    
    # plot linear model if significant
    if(sMantel$signif < 0.05){
      gg_list[[i]] <-
        ggplot(dist_mergeSUB, aes(.data[[metricDIST]], .data[[metricDIV]])) + #, color = Seychelles
        geom_smooth(method = "lm", color = "grey50", se = F) +
        geom_point(aes(shape = Seychelles), show.legend = F, color = color_y) +
        scale_shape_manual(values = c(19, 1)) +
        annotate('text',
                 x=min(dist_mergeSUB[[metricDIST]]), y=max(dist_mergeSUB[[metricDIV]]),
                 hjust = 0, vjust = 1,
                 label=paste0("r = ", sprintf("%.3f", sMantel$statistic), "\np = ", sMantel$signif)) +
        theme_light() +
        ylab(metricLAB) + xlab('Environmental distance') +
        theme(
          axis.title.x = element_text(colour = color_x),
        )
    }
    
    # no linear model if not significant
    if(sMantel$signif > 0.05 | (locations == "allsites" & metricDIV == "Fst")){
      gg_list[[i]] <-
        ggplot(dist_mergeSUB, aes(.data[[metricDIST]], .data[[metricDIV]], shape = Seychelles)) + #
        geom_point(aes(shape = Seychelles), show.legend = F, color = color_y) +
        scale_shape_manual(values = c(19, 1)) +
        annotate('text',
                 x=min(dist_mergeSUB[[metricDIST]]), y=max(dist_mergeSUB[[metricDIV]]),
                 hjust = 0, vjust = 1,
                 label=paste0("r = ", sprintf("%.3f", sMantel$statistic), "\np = ", sMantel$signif)) +
        theme_light() +
        ylab(metricLAB) + xlab('Environmental distance') +
        theme(
          axis.title.x = element_text(colour = color_x),
        )
    }
    
    IBE_table <- 
      IBE_table %>% 
      rbind(tibble(locations = locations,
                   metric = metricDIV,
                   isolation_by = metricDIST,
                   Mantel_r = sMantel$statistic,
                   Mantel_pval = sMantel$signif,
                   Procruste_r = sqrt(1-sProtest$ss),
                   Procruste_pval = sProtest$signif))
    
    i = i+1
  }
}


### {FIGURE S10} ####
gg_grob <-
  patchwork::wrap_plots(gg_list, ncol = 2, tag_level = "new") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
plot(gg_grob)

ggsave(gg_grob, width = 10, height = 15, dpi = 500,
       filename = paste0("results/4_continuity/_S10_isolation_by_", metricDIST, "_", comm_delin, "_2.png"))


### {TABLE S5} ####
IBE_table$Mantel_signif <- ifelse(IBE_table$Mantel_pval < 0.001, "***", 
                                  ifelse(IBE_table$Mantel_pval < 0.01, "**", 
                                         ifelse(IBE_table$Mantel_pval < 0.05, "*", "NS")))
IBE_table$Procruste_signif <- ifelse(IBE_table$Procruste_pval < 0.001, "***", 
                                     ifelse(IBE_table$Procruste_pval < 0.01, "**", 
                                            ifelse(IBE_table$Procruste_pval < 0.05, "*", "NS")))
IBE_table <-
  IBE_table %>%
  mutate(Mantel_r = paste(round(Mantel_r, 3), Mantel_signif)) %>%
  mutate(Procruste_r = paste(round(Procruste_r, 3), Procruste_signif)) %>%
  mutate(Mantel_r = gsub(" NS", "", Mantel_r)) %>%
  mutate(Procruste_r = gsub(" NS", "", Procruste_r)) %>%
  dplyr::select(-c(Mantel_pval, Procruste_pval, Mantel_signif, Procruste_signif)) %>%
  dplyr::arrange(locations)

write_csv(IBE_table, paste0("results/4_continuity/_S5_isolation_by_", metricDIST, "_", comm_delin, ".csv"))



# ---- 2. SGDC table ----
## table ----
library(boot)

metricGD = "Fst"

library(ecodist)        # MRM

SGDC_beta <- tibble()
i=1 # initialisation
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
      
      fm <- as.formula(paste(metricSDcomm, " ~ ", metricGD, sep=""))
      sMantel3 <- ecodist::mantel(dist_mat[[metricSDcomm]] ~ dist_mat[[metricGD]], nperm = 9999, nboot = 9999, cboot = 0.95)        
      
      ### Pearson ----
      sPearson <- cor.test(dist_merge[[metricSDcomm]], dist_merge[[metricGD]],
                           alternative = c("two.sided"), method = "pearson")
      sPearsonPERM <- jmuOutlier::perm.cor.test(dist_merge[[metricSDcomm]], dist_merge[[metricGD]], 
                                                alternative = c("two.sided"), method = "pearson", 
                                                num.sim = 9999)
      
      
      ### Procruste ----
      sProtest <- vegan::protest(dist_mat[[metricGD]], dist_mat[[metricSDcomm]], permutations = 9999)
      sProtest
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
      SGDC_beta[i, "Pearson_r"] <- sPearson$estimate
      SGDC_beta[i, "Pearson_pval"] <- sPearsonPERM$p.value
      SGDC_beta[i, "Procruste_r"] <- sqrt(1-sProtest$ss)
      SGDC_beta[i, "Procruste_pval"] <- sProtest$signif
      
      cat(comm_delin, locations, comm, metricSD, "\n")
      i=i+1
    }
  }
}


SGDC_beta$Mantel_signif <- ifelse(SGDC_beta$Mantel_pval < 0.001, "***", 
                                  ifelse(SGDC_beta$Mantel_pval < 0.01, "**", 
                                         ifelse(SGDC_beta$Mantel_pval < 0.05, "*", "NS")))
SGDC_beta$Pearson_signif <- ifelse(SGDC_beta$Pearson_pval < 0.001, "***",
                                   ifelse(SGDC_beta$Pearson_pval < 0.01, "**",
                                          ifelse(SGDC_beta$Pearson_pval < 0.05, "*", "NS")))
SGDC_beta$Procruste_signif <- ifelse(SGDC_beta$Procruste_pval < 0.001, "***", 
                                     ifelse(SGDC_beta$Procruste_pval < 0.01, "**", 
                                            ifelse(SGDC_beta$Procruste_pval < 0.05, "*", "NS")))

SGDC_beta %>%
  write_csv("results/4_continuity/beta_SGDCs_table_CIlm_Pearson_Procruste_taxonomy.csv")

temp <-
  SGDC_beta %>%
  dplyr::filter(metricSD != "beta.jne") %>%
  mutate(Mantel_r = paste(round(Mantel_r, 3), Mantel_signif)) %>%
  mutate(Procruste_r = paste(round(Procruste_r, 3), Procruste_signif)) %>%
  mutate(Mantel_r = gsub(" NS", "", Mantel_r)) %>%
  mutate(Procruste_r = gsub(" NS", "", Procruste_r)) %>%
  dplyr::select(locations, metricGD, metricSD, taxonomic_scale, Mantel_r, Procruste_r) %>%
  dplyr::arrange(locations, metricSD)

temp
temp %>%
  write_csv("results/4_continuity/_beta_SGDCs_table_Mantel_Procruste.csv")





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
  scale_shape_manual(values=c(19, 19, 1), )+
  scale_color_manual(values = c("#60CEACFF", "#395D9CFF"))+ # "0B0405FF", 
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
  dplyr::filter(community_delineation == "taxonomy") %>% 
  dplyr::filter(metricSD == "beta.jtu") %>% 
  dplyr::filter(locations == "without Seychelles")

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







# ---- 4. Covariance decomposition ----

## load ----
source("scripts/4_continuity/sgdcs_decomposition_Lamy.R")

metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"

## subset Seychelles ----
dist_merge <- mat.subset(dist_merge, "Seychelles")

## plot ----
comm = names_communities[1]
gglist <- list()

tempDF <- tibble()

for (comm in names_communities[1:2]){
  metricSD = paste(comm, "beta.jtu", sep = ".")
  
  temp <- SGDC.decomp(SD = dist_merge[[metricSD]],
                      GD = dist_merge[[metricGD]],
                      FACTOR = dist_merge[,c("environment","seadist")])
  temp <- 
    temp$Contribution %>% 
    as_tibble(rownames = NA) %>% 
    rownames_to_column("variable") %>% 
    mutate(community = comm)
  
  # rename variables
  temp$variable <- gsub('Var\\(', '', temp$variable)
  temp$variable <- gsub('Cov\\(', '', temp$variable)
  temp$variable <- gsub('\\)', '', temp$variable)
  temp$variable <- gsub('environment', 'Environment', temp$variable)
  temp$variable <- gsub('seadist', 'Distance', temp$variable)
  temp$variable <- gsub(',', ' ×', temp$variable)
  temp$variable <- gsub('SGDC', 'β-SGDC', temp$variable)
  
  
  temp$variable <- factor(temp$variable, levels = temp$variable[c(3,2,1,4,5)])
  
  
  ### {FIGURE 5} ####
  gglist[[comm]] <-
    ggplot(temp, aes(y = variable, x = Contribution, fill = variable)) +
    geom_bar(stat = "identity") +
    xlab('') + ylab('') +
    scale_fill_manual(values = c("#ffb734", "#9ECAE1", "#A1D99B", "grey", "grey20")) +
    theme_light() +
    theme(legend.position = "none",
          text = element_text(size = 15))
  
}

gg_grob <-
  patchwork::wrap_plots(gglist, nrow = 2, tag_level = "new") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
plot(gg_grob)

ggsave(paste0("results/4_continuity/_5_variance_decompositionPERSO_beta_SGDCs_noSeychelles_", 
              names_communities[1], "_", names_communities[2], "_", comm_delin, "_bdmean_r.png"),
       width = 8, height = 7, dpi = 500)
ggsave(paste0("results/4_continuity/_5_variance_decompositionPERSO_beta_SGDCs_noSeychelles_", 
              names_communities[1], "_", names_communities[2], "_", comm_delin, "_bdmean_r.pdf"),
       width = 8, height = 7)




# ---- 6. dbRDA geo + env ----
## dbRDA GD table ----
library(adespatial)
library(spdep)
data(oribatid)
nbtri <- tri2nb(as.matrix(coord_sites))
sc.tri <- scores.listw(nb2listw(dbMEM_seadist, style = "B"))
summary(sc.tri)
sc.tri


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
  
  ### a. seadist dbMEM ----
  dbMEM_seadist <- adespatial::dbmem(dist_mat$seadist, store.listw = TRUE, MEM.autocor = "positive") # MEM.autocor = "non-null",
  
  ## keep MEMs with Moran's I > 0.25
  # Moran's I values are proportional to Eigenvalues 
  MoranI_seadist <- moran.randtest(dbMEM_seadist, nrepet = 999)
  print(which(MoranI_seadist$obs > 0.25))
  
  dbMEM_seadist <- 
    dbMEM_seadist %>% 
    as.data.frame() %>% 
    dplyr::select(as.numeric(which(MoranI_seadist$obs > 0.25)))
  
  
  ### b. environment PCA ----
  ## scale
  envt_site_scaled <-
    envt_site %>%
    column_to_rownames("site") %>%
    scale(center = T, scale = T) %>%
    as.data.frame()
  
  ## PCA
  PCA_envt <- dudi.pca(envt_site_scaled, scannf = F, nf = 2)
  PCA_envtaxis <- as.data.frame(PCA_envt$li) # The principal components can be found in the $x matrix
  
  ## keep fist 2 axis, explains ~95% of variation
  PCA_envtaxis <- PCA_envtaxis[, 1:2]
  colnames(PCA_envtaxis) <- c("PC1", "PC2")
  
  
  for (metricGD in c("Fst", "GstPP.hed", "D.Jost")){ 
    
    ### c. best model table ----
    X = cbind(dbMEM_seadist, PCA_envtaxis)
    ade4::s.value(dfxy = coord_sites, z = X[, c("MEM1")])
    
    dbrda_GD0 <- vegan::capscale(dist_mat[[metricGD]] ~ 1, data = X) # null model, only intercept
    dbrda_GDall <- vegan::capscale(dist_mat[[metricGD]] ~ ., data = X)
    
    dbrda_GDsel <- ordistep(dbrda_GD0, scope = formula(dbrda_GDall), 
                            direction="both", trace = FALSE) 
    # dbrda_GDsel$anova
    anova_GDsel <- anova.cca(dbrda_GDsel, permutations = 9999)
    
    
    if(!is.null(dbrda_GDsel$anova)){
      dbRDAgd <-
        dbRDAgd %>%
        rbind(tibble(location = locations,
                     Model = as.character(dbrda_GDsel$call)[2],
                     R2 = RsquareAdj(dbrda_GDsel)$r.squared, 
                     adjR2 = RsquareAdj(dbrda_GDsel)$adj.r.squared,
                     F = anova_GDsel[1,3],
                     pval = anova_GDsel[1,4])) %>% 
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
    
    
    ### d. varpart ----
    ## with ALL explanatory variables
    vpGD <- varpart(dist_mat[[metricGD]], dbMEM_seadist, PCA_envtaxis)
    vpGD
    showvarparts(2, bg = c("#9ECAE1", "#A1D99B"))
    plot(vpGD, bg = c("#9ECAE1", "#A1D99B"), Xnames = c("Distance", "Environment"),
         main = metricGD)
    
    
    ## Varpart computation by hand : sligthly different from varpart plot...
    RDAabc <- vegan::capscale(dist_mat[[metricGD]] ~ ., data = cbind(dbMEM_seadist, PCA_envtaxis))
    RDAa <- vegan::capscale(dist_mat[[metricGD]] ~ MEM1 + Condition(PC1 + PC2), data = cbind(dbMEM_seadist, PCA_envtaxis))
    RDAb <- vegan::capscale(dist_mat[[metricGD]] ~ PC1 + PC2 + Condition(MEM1), data = cbind(dbMEM_seadist, PCA_envtaxis))
    
    a <- RsquareAdj(RDAa)$r.squared
    b <- RsquareAdj(RDAb)$r.squared
    c <- RsquareAdj(RDAabc)$r.squared - a - b
    
    varpartGD <- 
      varpartGD %>% 
      rbind(tibble(location = locations,
                   response_variable = metricGD,
                   explanatory_variable = c("Distance", "Environment", "Shared", "Residuals", "Total"),
                   variation_explained = c(a, b, c, 1-a-b-c, a+b+c),
                   pval = c(anova(RDAa)[1,4], anova(RDAb)[1,4], NA, NA, anova(RDAabc)[1,4]))
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

varpartGDwide <- 
  varpartGD %>% 
  mutate(variation_explained = paste(round(variation_explained, 4), signif)) %>% 
  mutate(variation_explained = gsub(" NA", "", variation_explained)) %>% 
  dplyr::select(-c(pval, signif)) %>% 
  pivot_wider(names_from = explanatory_variable, values_from = variation_explained)
varpartGDwide

# save
dbRDAgd %>% write_csv("results/4_continuity/dbRDA_models_beta_gd_ALL.csv")
varpartGDwide %>% write_csv("results/4_continuity/dbRDA_varpart_beta_gd_ALL.csv")





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
  
  ### a. dbMEM on seadist ----
  dbMEM_seadist <- adespatial::dbmem(dist_mat$seadist, store.listw = TRUE, MEM.autocor = "positive") # MEM.autocor = "non-null",
  
  ## keep MEMs with Moran's I > 0.25
  MoranI_seadist <- moran.randtest(dbMEM_seadist, nrepet = 999)
  print(which(MoranI_seadist$obs > 0.25))
  dbMEM_seadist <- 
    dbMEM_seadist %>% 
    as.data.frame() %>% 
    dplyr::select(as.numeric(which(MoranI_seadist$obs > 0.25)))
  
  
  ### b. environment PCA (X) ----
  ## scale
  envt_site_scaled <-
    envt_site %>% 
    column_to_rownames("site") %>% 
    scale(center = T, scale = T) %>% 
    as.data.frame()
  
  ## PCA
  PCA_envt <- dudi.pca(envt_site_scaled, scannf = F, nf = 2)
  PCA_envtaxis <- as.data.frame(PCA_envt$li) # The principal components can be found in the $x matrix
  
  ## keep fist 2 axis, explains ~95% of variation
  PCA_envtaxis <- PCA_envtaxis[, 1:2]
  colnames(PCA_envtaxis) <- c("PC1", "PC2")
  
  
  
  for (comm in names_communities){
    
    for (metricSD in c("beta.jac", "beta.jtu")){ #"beta.jne"
      
      metricSDcomm <- paste(comm, metricSD, sep = ".")
      
      ### c. best model table ----
      X = cbind(dbMEM_seadist, PCA_envtaxis)
      dbrda_SD0 <- vegan::capscale(dist_mat[[metricSDcomm]] ~ 1, data = X) # null model, only intercept
      dbrda_SDall <- vegan::capscale(dist_mat[[metricSDcomm]] ~ ., data = X)
      
      dbrda_SDsel <- ordistep(dbrda_SD0, scope = formula(dbrda_SDall), 
                              direction="both", trace = FALSE) 
      anova_SDsel <- anova.cca(dbrda_SDsel, permutations = 9999)
      
      if(!is.null(dbrda_SDsel$anova)){
        dbRDAsd <-
          dbRDAsd %>% 
          rbind(tibble(location = locations,
                       Model = as.character(dbrda_SDsel$call)[2],
                       R2 = RsquareAdj(dbrda_SDsel)$r.squared, 
                       adjR2 = RsquareAdj(dbrda_SDsel)$adj.r.squared,
                       F = anova_SDsel[1,3],
                       pval = anova_SDsel[1,4])) %>% 
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
      
      ## Varpart computatio by hand : sligthly different from varpart plot...
      RDAabc <- vegan::capscale(dist_mat[[metricSDcomm]] ~ ., data = cbind(dbMEM_seadist, PCA_envtaxis))
      RDAa <- vegan::capscale(dist_mat[[metricSDcomm]] ~ MEM1  + Condition(PC1 + PC2), data = cbind(dbMEM_seadist, PCA_envtaxis))
      RDAb <- vegan::capscale(dist_mat[[metricSDcomm]] ~ PC1 + PC2 + Condition(MEM1), data = cbind(dbMEM_seadist, PCA_envtaxis))
      
      a <- RsquareAdj(RDAa)$adj.r.squared
      b <- RsquareAdj(RDAb)$adj.r.squared
      c <- RsquareAdj(RDAabc)$adj.r.squared - a - b
      
      varpartSD <- 
        varpartSD %>% 
        rbind(tibble(location = locations,
                     response_variable = metricSDcomm,
                     explanatory_variable = c("Distance", "Environment", "Shared", "Residuals", "Total"),
                     variation_explained = c(a, b, c, 1-a-b-c, a+b+c),
                     pval = c(anova(RDAa)[1,4], anova(RDAb)[1,4], NA, NA, anova(RDAabc)[1,4]))
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
dbRDAsd %>% write_csv("results/4_continuity/dbRDA_models_beta_sd_ALL.csv")
varpartSDwide %>% write_csv("results/4_continuity/dbRDA_varpart_beta_sd_ALL.csv")



