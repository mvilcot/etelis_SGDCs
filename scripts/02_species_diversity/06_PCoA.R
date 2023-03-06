
library(glue)

## ---- setup beta data ----
# get distance matrix

level = "site"
sd_beta <- readRDS(paste0("results/02_species_diversity/sd_list_pairwise_", level, ".RDS"))



## ---- pcoa analysis ----

gg_list <- list()
i=1
for (metricSD in names(sd_beta[[1]])){
  for (comm in names(sd_beta)){
    
    # read distance matrix
    mat_SDbeta <- sd_beta[[comm]][[metricSD]]
    
    # run pcoa
    pcoa <- dudi.pco(quasieuclid(mat_SDbeta), scannf=FALSE, nf=8)
    
    # eigenvalues
    percent_explained <- pcoa$eig / sum(pcoa$eig) * 100
    
    # set labels
    pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
    labels <- c(glue("PCo Axis 1 ({pretty_pe[1]}%)"),
                glue("PCo Axis 2 ({pretty_pe[2]}%)"))
    
    # plot
    gg_list[[i]] <- 
      ggplot(pcoa$li, aes(x=A1, y=A2, color=rownames(pcoa$li))) +
      geom_point() +
      labs(x=labels[1], y=labels[2]) +
      geom_text(label = rownames(pcoa$li)) +
      ggtitle(comm, subtitle = metricSD) +
      theme(legend.position="none")
    
    i=i+1
  }
  
}

gg_grob <- arrangeGrob(grobs = gg_list, 
          ncol = length(names(sd_beta)), 
          nrow = length(names(sd_beta[[1]])))

ggsave(gg_grob, width = 20, height = 10, 
       filename = "results/02_species_diversity/PCoA_beta_species_diversity.png")




## ---- load data ----
level = "site"

gd_global <- read.csv(paste0("results/01_genetic_diversity/gd_table_global_", level, ".csv"))
sd_global <- read.csv(paste0("results/02_species_diversity/sd_table_global.csv"))

gd_alpha <- read.csv(paste0("results/01_genetic_diversity/gd_table_", level, ".csv"))
sd_alpha <- read.csv(paste0("results/02_species_diversity/sd_table_", level, ".csv"))


## ---- handle data ----
# remove sites
gd_alpha <- 
  gd_alpha %>% 
  # filter(site != "Hawaii") %>% 
  filter(site != "Seychelles")

sd_alpha <- 
  sd_alpha %>% 
  # filter(site != "Hawaii") %>% 
  filter(site != "Seychelles")

# scale metrics
gd_alpha <- 
  gd_alpha %>%
  mutate_at(c("Ho", "Hs"), scale)

sd_alpha <- 
  sd_alpha %>%
  group_by(community) %>%
  mutate_at("richness_site", scale)


# merge tables
table_alpha <- 
  gd_alpha %>% 
  left_join(sd_alpha, by = level) %>% 
  select(site, community, Hs, richness_site) %>% 
  pivot_longer(cols = c("Hs", "richness_site"), names_to = "metric")

ggplot(sd_alpha, aes(x=site, y=richness_site, color=community)) + 
  geom_boxplot() 