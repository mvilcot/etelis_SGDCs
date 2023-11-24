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
# comm_delin <- "taxonomy_PCADAPTallsites"


# list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))
names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria/misc", "Teleostei")


## load ----
dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))

names(dist_mat)


## subset sites ----
dist_mat <- lapply(dist_mat, function(x) mat.subset(x, "Seychelles"))
dist_merge <- mat.subset(dist_merge, "Seychelles")











# ---- 0. check assumptions ----
model <- cor(dist_merge[-c(1:3)])
corrplot::corrplot(model)

hist(dist_mat[[metricSD]])
hist(dist_mat[[metricGD]]) # GD not normal

shapiro.test(dist_mat[[metricSD]])
shapiro.test(dist_mat[[metricGD]]) # GD not normal
qqnorm(dist_mat[[metricSD]])
qqnorm(dist_mat[[metricGD]])


# ---- 1. IBD/IBE + SGDCs ----

# parameters
comm <- "Lutjanidae"
metricSD <- paste(comm, "beta.jtu", sep = ".")
metricGD <- "Fst"
metricDIST <- "seadist"

# # MRM
# sMRM_IBDsd <- MRM(dist_mat[[metricSD]] ~ dist_mat[[metricDIST]], nperm = 9999)
# sMRM_IBDgd <- MRM(dist_mat[[metricGD]] ~ dist_mat[[metricDIST]], nperm = 9999)
# sMRM_SGDC <- MRM(dist_mat[[metricSD]] ~ dist_mat[[metricGD]], nperm = 9999)
# # MRM(dist_merge[[metricSD]] ~ dist_merge[[metricDIST]], nperm = 9999) # same with distance matrix or long df

# Mantel
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
           label=paste0("r Mantel = ", sprintf("%.3f", sMantel_IBDsd$statistic), ", p = ", sMantel_IBDsd$signif)) +
  # "\nR² MRM = ", round(sMRM_IBDsd$r.squared[["R2"]], 4), ", p = ", round(sMRM_IBDsd$r.squared[["pval"]], 5))) +
  theme_light()

# genetic IBD
ggGD <-
  ggplot(dist_merge) +
  geom_point(aes(.data[[metricDIST]], .data[[metricGD]], color = Seychelles), show.legend = F) +
  scale_color_brewer(palette="Dark2") +
  annotate('text',
           x=min(dist_merge[[metricDIST]]), y=max(dist_merge[[metricGD]]),
           hjust = 0, vjust = 1,
           label=paste0("r Mantel = ", sprintf("%.3f", sMantel_IBDgd$statistic), ", p = ", sMantel_IBDgd$signif)) +
  # "\nr MRM = ", round(sMRM_IBDgd$r.squared[["R2"]], 4), ", p = ", round(sMRM_IBDgd$r.squared[["pval"]], 5))) +
  theme_light()


# SGDCs
ggSGDCs <-
  ggplot(dist_merge) +
  geom_point(aes(.data[[metricSD]], .data[[metricGD]], color = Seychelles)) +
  scale_color_brewer(palette="Dark2") +
  annotate('text',
           x=min(dist_merge[[metricSD]]), y=max(dist_merge[[metricGD]]),
           hjust = 0, vjust = 1,
           label=paste0("r Mantel = ", sprintf("%.3f", sMantel_SGDC$statistic), ", p = ", sMantel_SGDC$signif)) +
  # "\nr MRM = ", round(sMRM_SGDC$r.squared[["R2"]], 4), ", p = ", round(sMRM_SGDC$r.squared[["pval"]], 5))) +
  theme_light()


# merge plots
ggSD + ggGD + ggSGDCs + #+ plot_annotation(title = comm)
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave(width = 14, height = 4,
       filename = paste0("results/4_continuity/IBD_beta_SGDC_noSeychelles_", metricSD, "_", metricGD, "_", metricDIST, "_", comm_delin, ".png"))



# ---- 2. SGDC table ----
## table ----
library(boot)
names_communities2 <- gsub("Eupercaria/misc", "Eupercaria_misc", names_communities)

metricGD = "Fst"

# create empty table
SGDC_beta <- tibble(community_delineation = NA,
                    locations = NA,
                    taxonomic_scale = NA,
                    metricGD = NA,
                    metricSD = NA,
                    r = NA,
                    rlwr = NA,
                    rupr = NA,
                    pval = NA,
)

# coeff_function <- function(formula, data, indices) {
#   d <- data[indices,] #allows boot to select sample
#   fit <- lm(formula, data=d) #fit regression model
#   return(summary(fit)$coefficients[2,1]) #return R-squared of model
# }
library(ecodist)        # MRM

i=1 # initialisation
for(comm_delin in comm_delin_list[1]) { #[1]

  for (locations in c("all sites", "without Seychelles")){

    dist_mat <- readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".RDS"))
    dist_merge <- read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))
    
    ## >>>> Eupercaria!!!!!!!! ----
    colnames(dist_merge) <- gsub("Eupercaria/misc", "Eupercaria_misc", colnames(dist_merge))
    names(dist_mat) <- gsub("Eupercaria/misc", "Eupercaria_misc", names(dist_mat))
    
    ## subset sites --
    if(locations == "without Seychelles"){
      dist_mat <- lapply(dist_mat, function(x){mat.subset(x, "Seychelles")})
      dist_merge <- mat.subset(dist_merge, "Seychelles")
    }

    for (comm in names_communities2){

      for (metric in c("beta.jac", "beta.jtu", "beta.jne")){

        metricSDcomm <- paste(comm, metric, sep = ".")

        stat <- vegan::mantel(dist_mat[[metricSDcomm]], dist_mat[[metricGD]], permutations = 9999)
        
        ## Get bootstraps on R²
        # fm <- as.formula(paste("scale(", metricSDcomm, ") ~ scale(", metricGD, ")", sep=""))
        # stat2 <- summary(lm(fm, data = dist_merge))
        # boots <- boot(data=dist_merge, statistic=coeff_function, R=1000, formula=fm)
        
        fm <- as.formula(paste(metricSDcomm, " ~ ", metricGD, sep=""))
        stat3 <- ecodist::mantel(dist_mat[[metricSDcomm]] ~ dist_mat[[metricGD]], nperm = 9999, nboot = 9999, cboot = 0.95)        
        # boots_ci <- boot.ci(boots, type="bca")
        # plot(reps)

        SGDC_beta[i,]$community_delineation <- comm_delin
        SGDC_beta[i,]$locations <- locations
        SGDC_beta[i,]$taxonomic_scale <- comm
        SGDC_beta[i,]$metricGD <- metricGD
        SGDC_beta[i,]$metricSD <- metric
        SGDC_beta[i,]$r <- stat3[1]# stat$statistic
        SGDC_beta[i,]$rlwr <- stat3[5] #stat2$coefficients[2,1] - 1.96*stat2$coefficients[2,2] #boots_ci$bca[4] 
        SGDC_beta[i,]$rupr <- stat3[6] # stat2$coefficients[2,1] + 1.96*stat2$coefficients[2,2] #boots_ci$bca[5]
        SGDC_beta[i,]$pval <- stat3[2] # stat$signif Ha: r>0


        cat(i, "\n")
        i=i+1
      }
    }
  }
}


SGDC_beta$signif <- ifelse(SGDC_beta$pval < 0.001, "***", 
                            ifelse(SGDC_beta$pval < 0.01, "**", 
                                   ifelse(SGDC_beta$pval < 0.05, "*", "NS")))

SGDC_beta %>%
  write_csv("results/4_continuity/beta_SGDCs_table_CImantel_taxonomy.csv")


## plot one community delin ----
SGDC_betasub <-
  SGDC_beta %>% 
  dplyr::filter(community_delineation == "taxonomy") %>% 
  dplyr::filter(locations == "without Seychelles") %>% 
  dplyr::filter(metricSD != "beta.jne")

SGDC_betasub$taxonomic_scale <- 
  factor(SGDC_betasub$taxonomic_scale,
         level = names_communities2)

SGDC_betasub$metricSD <- 
  factor(SGDC_betasub$metricSD,
         level = c("beta.jac", "beta.jtu")) #"beta.jne", 

ggplot(SGDC_betasub, aes(
                       x=r,
                       y=taxonomic_scale,
                       # shape=signif,
                       group=metricSD,
                       color=metricSD)) +
  geom_vline(xintercept = 0, linetype=2, color="grey20") +
  # geom_point(size = 2.5) +
  geom_pointrange(aes(x=r,
                      xmin=rlwr,
                      xmax=rupr),
                    size = 0.3, linewidth = 0.3,
                    position=position_dodge(width = 0.4)) +
  # scale_shape_manual(values=c(19, 19, 1), )+
  scale_color_manual(values = c("#60CEACFF", "#395D9CFF"))+ # "0B0405FF", 
  # scale_color_viridis(option="mako", discrete=T, end = 0.8) +
  xlim(c(-0.6, 0.9)) +
  ylab("") +
  xlab("Mantel r") +
  theme_light() +
  theme(legend.title=element_blank(), legend.position="top") +
  guides(shape = "none") +
  scale_y_discrete(limits=rev)

ggsave(paste0("results/4_continuity/beta_SGDCs_DWplot_taxonomy_noSeychelles_CImantel.png"),
       height = 4, width = 5, units = 'in', dpi = 300)



## plot with all comm delin ----
SGDC_beta$taxonomic_scale <- factor(SGDC_beta$taxonomic_scale,
                                    level = names_communities)

ggplot(SGDC_beta, aes(x=community_delineation,
                      y=r,
                      shape=signif,
                      group=taxonomic_scale,
                      color=taxonomic_scale)) +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(size = 2.5, position=position_dodge(width = 0.4)) +
  scale_shape_manual(values=c(19, 19, 19, 1))+
  scale_color_manual(values = c("#92e101", "#78c556", "#5fa8aa", "#458cff"))+
  facet_grid(vars(metricSD), vars(locations), scale="fixed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(paste0("results/4_continuity/beta_SGDCs_DWplot.png"),
       height = 8, width = 15, units = 'in', dpi = 300)



# ---- 3. SGDCs significant plot ----


SGDC_beta <- read_csv("results/4_continuity/beta_SGDCs_table.csv")


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
  
  dist_merge <-
    read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))
  
  if(SGDC_betaSUB[i,"locations"] == "without Seychelles"){
    dist_merge <- mat.subset(dist_merge, "Seychelles")
  }

  if(SGDC_betaSUB[i, "signif"] != "NS"){
    gg_list[[i]] <-
      ggplot(dist_merge, aes(.data[[metricSDcomm]], .data[[metricGD]])) +
      geom_smooth(method = "lm", color = "grey50", fill = "grey80") +
      geom_point(shape = 19, alpha = 0.8) +
      xlab(paste0(metricSD, " (", comm, ")")) +
      annotate('text',
               x = min(dist_merge[[metricSDcomm]]), y = 1.1 * max(dist_merge[[metricGD]]),
               hjust = 0, vjust = 1,
               label = paste0("r Mantel = ", round(SGDC_betaSUB[i, "r"], 4), "\np = ", SGDC_betaSUB[i, "pval"])) +
      theme_light()
  }
  if(SGDC_betaSUB[i, "signif"] == "NS"){
    gg_list[[i]] <-
      ggplot(dist_merge, aes(.data[[metricSDcomm]], .data[[metricGD]])) +
      geom_point(shape = 19, alpha = 0.8) +
      xlab(paste0(metricSD, " (", comm, ")")) +
      annotate('text',
               x = min(dist_merge[[metricSDcomm]]), y = 1.1 * max(dist_merge[[metricGD]]),
               hjust = 0, vjust = 1,
               label = paste0("r Mantel = ", round(SGDC_betaSUB[i, "r"], 4), "\np = ", SGDC_betaSUB[i, "pval"])) +
      theme_light()
  }
  
}

gg_grob <-
  patchwork::wrap_plots(gg_list, tag_level = "new") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
plot(gg_grob)
ggsave(gg_grob, width = 10, height = 9,
       filename = paste0("results/4_continuity/beta_SGDCs_plot_significant.png"))








# ---- 4. variance decomposition----

## load ----
source("scripts/4_continuity/sgdcs_decomposition_Lamy.R")
## >>>> à adapter à du beta... avec phytools::multi.mantel

comm_delin <- comm_delin_list[1]

dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))

comm = names_communities[2]
metricSD = paste(comm, "beta.jtu", sep = ".")
metricGD = "Fst"

## subset Seychelles ----
dist_merge <- mat.subset(dist_merge, "Seychelles")

## plot ----
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



# perso plot
temp <- SGDC.decomp(SD = dist_merge[[metricSD]],
                    GD = dist_merge[[metricGD]],
                    FACTOR = dist_merge[,c("environment","seadist")])
temp <- 
  temp$Contribution %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("variable")
temp$variable <- factor(temp$variable, levels = temp$variable)
ggplot(temp, aes(y = variable, x = Contribution, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#A1D99B", "#9ECAE1", "#ffb734", "grey", "grey20")) +
  theme_light() +
  theme(legend.position = "none")
ggsave(paste0("results/4_continuity/variance_decompositionPERSO_beta_SGDCs_noSeychelles_", comm, "_", comm_delin, "_bdmean.png"),
       width = 10, height = 5)
    
# ---- 5. MRM ----

## MRM gd table ----
comm_delin <- comm_delin_list[1]

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


### GD + SD ----
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

ggsave(paste0("results/4_continuity/DWplotPERSO_MRM_beta_combined_envt_seadist.png"),
         height = 8, width = 8, units = 'in', dpi = 500)




# ### GD ----
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
# ### SD ----
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
      xlab(paste0("residuals(", metricSD, ") (", comm, ")")) +
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
      xlab(paste0("residuals(", metricSD, ") (", comm, ")")) +
      ylab(paste0("residuals(", metricGD, ")")) +
      annotate('text',
               x = min(df$SDres), y = 1.5*max(df$GDres),
               hjust = 0, vjust = 1,
               label = paste0("r Mantel = ", round(stats_res$statistic, 4), "\np = ", stats_res$signif)) +
      theme_light()
  }
}

gg_grob <-
  patchwork::wrap_plots(gglist, tag_level = "new") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
plot(gg_grob)

ggsave(width = 10, height = 9,
       filename = paste0("results/4_continuity/beta_SGDCs_plot_significant_RESIDUALS.png"))





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


