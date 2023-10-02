
# ---- parameters ----
## communtity delineation
comm_delin = "taxonomy"
# comm_delin = "taxonomy_depth1_crosses45-400m"
# comm_delin = "taxonomy_depth2_contains45-400m"
# comm_delin = "taxonomy_depth3_within45-400m"
# comm_delin = "taxonomy_env_reef-associated"


list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

# parameters
comm = names(list_communities)[2]
comm
metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"
metricDIST = "seadist"


# ---- load ----
dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))

names(dist_mat)


# ---- subset sites ----
loc = "Seychelles"

for(i in 1:length(dist_mat)){
  dist_mat[[i]] <- as.matrix(dist_mat[[i]])
  dist_mat[[i]] <- dist_mat[[i]][!grepl(loc, rownames(dist_mat[[i]])),
                                 !grepl(loc, colnames(dist_mat[[i]]))]
  dist_mat[[i]] <- as.dist(dist_mat[[i]])
}

dist_merge <- dist_merge[!grepl(loc, dist_merge$site),]


# 
# # ---- 0. check assumptions ----
# model <- cor(dist_merge[-c(1:3)])
# corrplot::corrplot(model)
# 
# hist(dist_mat[[metricSD]])
# hist(dist_mat[[metricGD]]) # GD not normal
# 
# shapiro.test(dist_mat[[metricSD]])
# shapiro.test(dist_mat[[metricGD]]) # GD not normal
# 
# 

# ---- 1. IBD + SGDCs ----

# # MRM
# sMRM_IBDsd <- MRM(dist_mat[[metricSD]] ~ dist_mat[[metricDIST]], nperm = 9999)
# sMRM_IBDgd <- MRM(dist_mat[[metricGD]] ~ dist_mat[[metricDIST]], nperm = 9999)
# sMRM_SGDC <- MRM(dist_mat[[metricSD]] ~ dist_mat[[metricGD]], nperm = 9999)
# # MRM(dist_merge[[metricSD]] ~ dist_merge[[metricDIST]], nperm = 9999) # same with distance matrix or long df

# Mantel
#### !!!!!!!!!!! REMOVE values=0...... !!! ##################
sMantel_IBDsd <- vegan::mantel(dist_mat[[metricSD]], dist_mat[[metricDIST]], permutations = 9999)
sMantel_IBDgd <- vegan::mantel(dist_mat[[metricGD]], dist_mat[[metricDIST]], permutations = 9999)
sMantel_SGDC <- vegan::mantel(dist_mat[[metricSD]], dist_mat[[metricGD]], permutations = 9999)

# # distance decay models
# sDecay_IBDsd <- decay.model(dist_mat[[metricSD]], dist_mat[[metricDIST]], y.type="dissim", model.type="pow", perm=999)
# sDecay_IBDgd <- decay.model(dist_mat[[metricGD]], dist_mat[[metricDIST]], y.type="dissim", model.type="pow", perm=999)
#
# plot(dist_mat[[metricDIST]], dist_mat[[metricSD]], pch = 16)
# plot.decay(sDecay_IBDsd, col="green", remove.dots=T, add=T, lwd=2)
#
# plot(dist_mat[[metricDIST]], dist_mat[[metricGD]], pch = 16)
# plot.decay(sDecay_IBDgd, col="green", remove.dots=T, add=T, lwd=2)



# ---- plot ---- #
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
           label=paste0("r Mantel = ", round(sMantel_IBDsd$statistic, 4), ", p = ", round(sMantel_IBDsd$signif, 5))) +
  # "\nR² MRM = ", round(sMRM_IBDsd$r.squared[["R2"]], 4), ", p = ", round(sMRM_IBDsd$r.squared[["pval"]], 5))) +
  labs(title = "Species IBD")

# genetic IBD
ggGD <-
  ggplot(dist_merge) +
  geom_point(aes(.data[[metricDIST]], .data[[metricGD]], color = Seychelles), show.legend = F) +
  scale_color_brewer(palette="Dark2") +
  annotate('text',
           x=min(dist_merge[[metricDIST]]), y=max(dist_merge[[metricGD]]),
           hjust = 0, vjust = 1,
           label=paste0("r Mantel = ", round(sMantel_IBDgd$statistic, 4), ", p = ", round(sMantel_IBDgd$signif, 5))) +
  # "\nr MRM = ", round(sMRM_IBDgd$r.squared[["R2"]], 4), ", p = ", round(sMRM_IBDgd$r.squared[["pval"]], 5))) +
  labs(title = "Genetics IBD")


# SGDCs
ggSGDCs <-
  ggplot(dist_merge) +
  geom_point(aes(.data[[metricSD]], .data[[metricGD]], color = Seychelles)) +
  scale_color_brewer(palette="Dark2") +
  annotate('text',
           x=min(dist_merge[[metricSD]]), y=max(dist_merge[[metricGD]]),
           hjust = 0, vjust = 1,
           label=paste0("r Mantel = ", round(sMantel_SGDC$statistic, 4), ", p = ", round(sMantel_SGDC$signif, 4))) +
  # "\nr MRM = ", round(sMRM_SGDC$r.squared[["R2"]], 4), ", p = ", round(sMRM_SGDC$r.squared[["pval"]], 5))) +
  labs(title = "SGDC")


# merge plots
ggSD + ggGD + ggSGDCs + plot_annotation(title = comm)
ggsave(width = 20, height = 6,
       filename = paste0("results/4_continuity/IBD_beta_SGDC_noSeychelles_", metricSD, "_", metricGD, "_", metricDIST, "_", comm_delin, ".png"))



# ---- 2. SGDCs by community ----

gg_list <- list()
stat_Mantel <- 
  stat_MRM <- 
  stat_LM <- 
  list()
metricSD = "beta.jtu"

for (comm in names(list_communities)){
  metricSDcomm = paste(comm, metricSD, sep = ".")
  print(metricSDcomm)
  
  
  ## statistics ----
  # Mantel test
  stat_Mantel[[comm]] <- vegan::mantel(dist_mat[[metricSDcomm]], dist_mat[[metricGD]], permutations = 9999)
  
  # MRM
  stat_MRM[[comm]] <- MRM(dist_mat[[metricSDcomm]] ~ dist_mat[[metricGD]], nperm = 9999)
  
  # LM
  stat_LM[[comm]] <- lm(dist_merge[[metricSDcomm]] ~ dist_merge[[metricGD]])
  
  
  ## plot ----
  gg_list[[comm]] <-
    ggplot(dist_merge, aes(.data[[metricSDcomm]], .data[[metricGD]])) +
    geom_point() +
    xlab(metricSDcomm) +
    # geom_text_repel(aes(label = sites), size=3) +
    # ggtitle(paste0(patt)) +
    annotate('text',
             x=min(dist_merge[[metricSDcomm]]), y=max(dist_merge[[metricGD]]),
             hjust = 0, vjust = 1,
             label=paste0("r Mantel = ", round(stat_Mantel[[comm]]$statistic, 4), ", p = ", round(stat_Mantel[[comm]]$signif, 4))) +
    # "\nr MRM = ", round(stat_MRM[[comm]]$r.squared[["R2"]], 4), ", p = ", round(stat_MRM[[comm]]$r.squared[["pval"]], 5))) +
    labs(title = comm)
  
  
  # ggsave(gg, device = "png", height = 4, width = 6, units = "in",
  #        filename = paste0("results/04_continuity/SGDCs_beta_", patt, "_", level, "_", metricGD, "_", metricSD, "_", comm, ".png"))

  
}

gg_grob <- arrangeGrob(grobs = gg_list, ncol=4)
plot(gg_grob)
ggsave(gg_grob, width = 20, height = 5,
       filename = paste0("results/4_continuity/beta_SGDCs_noSeychelles_", level, "_", metricGD, "_", metricSD, "_",  comm_delin, ".png"))


# ---- 3. SGDC table ----

i=1
metricGD = "Fst"
SGDC_table <- tibble(community_delineation = NA,
                     locations = NA,
                     taxonomic_scale = NA,
                     metricGD = NA,
                     metricSD = NA,
                     r = NA,
                     pval = NA
                     )
for(comm_delin in comm_delin_list) {
  for (locations in c("allsites", "noSeychelles")){
    
    dist_mat <-
      readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".RDS"))
    
    names(dist_mat)
    
    
    ## subset sites --
    if(locations == "noSeychelles"){
      for(j in 1:length(dist_mat)){
        dist_mat[[j]] <- as.matrix(dist_mat[[j]])
        dist_mat[[j]] <- dist_mat[[j]][!grepl("Seychelles", rownames(dist_mat[[j]])),
                                       !grepl("Seychelles", colnames(dist_mat[[j]]))]
        dist_mat[[j]] <- as.dist(dist_mat[[j]])
      }
    }
    
    for (comm in names(list_communities)){
      
      for (metric in c("beta.jac", "beta.jtu", "beta.jne")){
        metricSD <- paste(comm, metric, sep = ".")

        stat <- vegan::mantel(dist_mat[[metricSD]], dist_mat[[metricGD]], permutations = 9999)
        
        SGDC_table[i,]$community_delineation <- comm_delin
        SGDC_table[i,]$locations <- locations
        SGDC_table[i,]$taxonomic_scale <- comm
        SGDC_table[i,]$metricGD <- metricGD
        SGDC_table[i,]$metricSD <- metric
        SGDC_table[i,]$r <- stat$statistic
        SGDC_table[i,]$pval <- stat$signif
        
        i=i+1
      }
    }
  }
}

SGDC_table %>% 
  write_csv("results/4_continuity/beta_SGDCs_table.csv")


## plot --
SGDC_table$taxonomic_scale <- factor(SGDC_table$taxonomic_scale, 
                                     level = unique(SGDC_table$taxonomic_scale))

SGDC_table$signif <- ifelse(SGDC_table$pval < 0.05, "*", "")

ggplot(SGDC_table, aes(x=community_delineation, 
                       y=r, 
                       shape=signif,
                       group=taxonomic_scale,
                       color=taxonomic_scale)) + 
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(size = 2.5, position=position_dodge(width = 0.4)) +
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values = c("#92e101", "#78c556", "#5fa8aa", "#458cff"))+
  facet_grid(vars(metricSD), vars(locations), scale="fixed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(paste0("results/4_continuity/beta_SGDCs_plot.png"),
       height = 8, width = 15, units = 'in', dpi = 300)





# ---- 4. variance decomposition----

source("scripts/4_continuity/sgdcs_decomposition_Lamy.R")

comm = names(list_communities)[2]
metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"


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



# ---- 5. MRM ----

## MRM gd table ----

MRMgd <- tibble::tibble()
for (metricGD in c("Fst", "GstPP.hed", "D.Jost")){
  for (locations in c("allsites", "noSeychelles")){
    dist_merge <-
      read_csv(paste0("results/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_taxonomy.csv"))
    
    if(locations == "noSeychelles"){dist_merge <- dist_merge[!grepl("Seychelles", dist_merge$site),]}
    
    statLM <- summary(lm(scale(dist_merge[[metricGD]]) ~
                           scale(dist_merge$environment) +
                           scale(dist_merge$seadist)))
    
    statMRM <- MRM(scale(dist_merge[[metricGD]]) ~
                     scale(dist_merge$environment) +
                     scale(dist_merge$seadist),
                   nperm = 9999)  
    
    tempLM <-
      data.frame(statLM$coefficients[-1,1:2]) %>%
      rownames_to_column("explanatory_variable")
    
    tempMRM <-
      data.frame(pval = statMRM$coef[-1,2]) %>%
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


MRMgd$explanatory_variable <- gsub('scale(dist_merge$', '', MRMgd$explanatory_variable, fixed = T)
MRMgd$explanatory_variable <- gsub(')', '', MRMgd$explanatory_variable, fixed = T)

MRMgd %>%
  write_csv("results/4_continuity/MRM_beta_gd_envt_seadist.csv")


## MRM sd table ----

comm_delin_list <- 
  c("taxonomy",
    "taxonomy_depth1_crosses45-400m",
    "taxonomy_depth2_contains45-400m",
    "taxonomy_depth3_within45-400m",
    "taxonomy_env_reef-associated")

MRMsd <- tibble::tibble()
for(comm_delin in comm_delin_list) {
  for (locations in c("allsites", "noSeychelles")){
    dist_merge <-
      read_csv(paste0("results/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".csv"))
    if(locations == "noSeychelles"){dist_merge <- dist_merge[!grepl("Seychelles", dist_merge$site),]}
    for (comm in names(list_communities)){
      
      for (metric in c("beta.jac", "beta.jtu", "beta.jne")){
        metricSD <- paste(comm, metric, sep = ".")
        statLM <- summary(lm(scale(dist_merge[[metricSD]]) ~
                               scale(dist_merge$environment) +
                               scale(dist_merge$seadist)))
        
        statMRM <- MRM(scale(dist_merge[[metricSD]]) ~
                         scale(dist_merge$environment) +
                         scale(dist_merge$seadist),
                       nperm = 9999)  
        
        tempLM <-
          data.frame(statLM$coefficients[-1,1:2]) %>%
          rownames_to_column("explanatory_variable")
        
        tempMRM <-
          data.frame(pval = statMRM$coef[-1,2]) %>%
          rownames_to_column("explanatory_variable")
        
        temp <-
          full_join(tempLM, tempMRM, by = "explanatory_variable") %>% 
          cbind(community_delineation = comm_delin,
                locations = locations,
                taxonomic_scale = comm,
                response_variable = metric,
                .)
        
        MRMsd <-
          MRMsd %>%
          rbind(temp) %>% 
          as_tibble()
      }
    }
  }
}


MRMsd$explanatory_variable <- gsub('scale(dist_merge$', '', MRMsd$explanatory_variable, fixed = T)
MRMsd$explanatory_variable <- gsub(')', '', MRMsd$explanatory_variable, fixed = T)

MRMsd %>%
  write_csv("results/4_continuity/MRM_beta_sd_envt_seadist.csv")


## cor explanatory variables ----
# model <- cor(dist_merge[c("geodist",
#                "seadist",
#                "leastcost",
#                "environment")])
# corrplot::corrplot(model)

cor(dist_merge[["seadist"]], dist_merge[["environment"]])




## DW plot by hand ----

### load
MRMgd <- read_csv("results/4_continuity/MRM_beta_gd_envt_seadist.csv")
MRMsd <- read_csv("results/4_continuity/MRM_beta_sd_envt_seadist.csv")

MRMgd$response_variable <- factor(MRMgd$response_variable, levels = unique(MRMgd$response_variable))
MRMsd$response_variable <- factor(MRMsd$response_variable, levels = unique(MRMsd$response_variable))

### GD
ggplot(MRMgd, aes(x=explanatory_variable, y=Estimate, group=response_variable, color=response_variable)) + 
  geom_hline(yintercept = 0,linetype=2) +
  geom_point(position=position_dodge(width = -0.4)) +
  geom_pointrange(aes(ymin=Estimate - 1.96*Std..Error,
                    ymax=Estimate + 1.96*Std..Error),
                position=position_dodge(width = -0.4)) +
  facet_wrap(~locations, scale="free") +
  ylim(c(-1,1)) +
  coord_flip()

ggsave(paste0("results/4_continuity/DWplotPERSO_MRM_beta_gd_envt_seadist.png"),
       height = 5, width = 12, units = 'in', dpi = 300)


### SD
# all combined
ggplot(MRMsd, aes(x=explanatory_variable, 
                  y=Estimate, 
                  shape=community_delineation, 
                  group=interaction(community_delineation, response_variable),
                  color=response_variable)) + 
  geom_hline(yintercept = 0,linetype=2) +
  geom_point(position=position_dodge(width = -1)) +
  geom_pointrange(aes(ymin=Estimate - 1.96*Std..Error,
                    ymax=Estimate + 1.96*Std..Error),
                position=position_dodge(width = -1)) +
  facet_grid(vars(taxonomic_scale), vars(locations), scale="free") +
  coord_flip()

ggsave(paste0("results/4_continuity/DWplotPERSO_MRM_beta_sd_envt_seadist_", comm_delin, ".png"),
       height = 8, width = 15, units = 'in', dpi = 300)


# beta.jtu separated

for(metric in unique(MRMsd$response_variable)){
  MRMsub <- 
    MRMsd %>% 
    filter(response_variable == metric)

  ggplot(MRMsub, aes(x=explanatory_variable, 
                    y=Estimate, 
                    shape=community_delineation,
                    # group=interaction(community_delineation, response_variable),
                    color=community_delineation)) + 
    geom_hline(yintercept = 0,linetype=2) +
    geom_point(position=position_dodge(width = -0.5)) +
    geom_pointrange(aes(ymin=Estimate - 1.96*Std..Error,
                        ymax=Estimate + 1.96*Std..Error),
                    position=position_dodge(width = -0.5)) +
    facet_grid(vars(taxonomic_scale), vars(locations), scale="free") +
    scale_x_discrete(limits = rev(levels(MRMsub$community_delineation))) +
    scale_color_viridis(option="mako", discrete=T)+
    coord_flip() +
    ggtitle(metric)
  
  ggsave(paste0("results/4_continuity/DWplotPERSO_MRM_beta_sd_envt_seadist_", comm_delin, "_", metric, ".png"),
         height = 8, width = 15, units = 'in', dpi = 300)
  
}





## DW plot automatic ----

### GD
models_gd <- list()
for (metricGD in c("Fst", "GstPP.hed", "D.Jost")){
  models_gd[[metricGD]] <- lm(scale(dist_merge[[metricGD]]) ~
                                scale(dist_merge$environment) +
                                scale(dist_merge$seadist))
  
}
dwplot(models_gd,
       vline = geom_vline(xintercept = 0,
                          colour = "grey60",
                          linetype = 2))

ggsave("results/4_continuity/DWplot_MRM_beta_gd_envt_seadist_noSeychelles.png",
       height = 5, width = 8, units = 'in', dpi = 300)

### SD
models_sd <- list()
for (metric in c("beta.jac", "beta.jtu", "beta.jne")){
  for (comm in names(list_communities)){
    
    metricSD <- paste(comm, metric, sep = ".")
    
    models_sd[[metricSD]] <- lm(scale(dist_merge[[metricSD]]) ~
                                  scale(dist_merge$environment) +
                                  scale(dist_merge$seadist))
  }
}
dwplot(models_sd,
       vline = geom_vline(xintercept = 0,
                          colour = "grey60",
                          linetype = 2))
ggsave(paste0("results/4_continuity/DWplot_MRM_beta_sd_envt_seadist_noSeychelles_", comm_delin, ".png"),
       height = 5, width = 8, units = 'in', dpi = 300)


# 6SGDCs residuals ----
# multi.mantel: same values as MRM ! but allows to take residuals

comm_delin = comm_delin_list[3]
list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

comm = names(list_communities)[2]
metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"

dist_mat <- readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".RDS"))

for(i in 1:length(dist_mat)){
  dist_mat[[i]] <- as.matrix(dist_mat[[i]])
  dist_mat[[i]] <- dist_mat[[i]][!grepl("Seychelles", rownames(dist_mat[[i]])),
                                 !grepl("Seychelles", colnames(dist_mat[[i]]))]
  dist_mat[[i]] <- as.dist(dist_mat[[i]])
}


## MRM ----
modelGD <-
  phytools::multi.mantel((dist_mat[[metricGD]]),
                         list(dist_mat$environment,
                              dist_mat$seadist),
                         nperm = 1000)

modelSD <-
  phytools::multi.mantel((dist_mat[[metricSD]]),
                         list(dist_mat$environment,
                              dist_mat$seadist),
                         nperm = 1000)

## residual correlation ----
mantel(modelGD$residuals, modelSD$residuals)








# ---- *** DRAFTS ----
## ---- *** SGDCs by distance phylo ----
# 
# SGDCs <- data.frame(community = names(list_communities))
# 
# 
# temp <- lapply(SGDCs$community, FUN = function(x){stat_Mantel[[x]]$statistic})
# names(temp) <- names(list_communities)
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
# 
