
# ---- Parameters ----
## communities delineation
comm_delin = "taxonomy"

list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))
names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")

# parameters
comm = names(list_communities)[2]
metricSD = "richness_site"


# ---- Load ----
level = "site"

sd_alpha <- read_csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))
sd_global <- read_csv(paste0("results/2_species_diversity/sd_table_global_", comm_delin, ".csv"))
sd_beta <- readRDS(paste0("results/2_species_diversity/sd_list_pairwise_site_", comm_delin, ".RDS"))

dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))


# ---- Summary ----
temp <-
  dist_merge %>% 
  dplyr::select(contains(c("beta.jtu", "beta.jac", "beta.jne")))

summary(temp)



# Violin plots by community ----
## alpha ----
gga <-
  sd_alpha %>% 
  mutate(site = factor(site, levels = levels(data_sites$site))) %>% 
  mutate(community = gsub("/misc", "", community)) %>% 
  mutate(community = factor(community, levels = names_communities)) %>% 
  ggplot(aes(community, richness_site)) +
  geom_violin(fill = "grey50", color = "grey50", scale = "width") +
  geom_boxplot(width=0.1) +
  xlab("") + ylab('α-diversity (species richness by site)') +
  theme_light()



## beta ----
ggb <-
  dist_merge %>% 
  dplyr::select(-c(Fst, GstPP.hed, D.Jost, jtu, jac, jne)) %>% 
  pivot_longer(cols = -c(site, site1, site2, geodist, seadist, environment), 
               names_to = "variable") %>% 
  separate(variable, c("community", "drop", "variable"), sep = "\\.") %>% 
  mutate(community = factor(community, levels = names_communities)) %>% 
  mutate(variable = paste0("β-", variable)) %>% 
  dplyr::filter(variable != "β-jne") %>% 
  ggplot(aes(community, value)) +
  geom_violin(aes(fill = variable), color = NA, scale = "width", position = position_dodge(0.8)) +
  geom_boxplot(aes(color = variable), width=0.1, position = position_dodge(0.8)) +
  scale_color_manual(values = c("grey20", "grey")) +
  scale_fill_manual(values = c("grey", "grey20")) +
  theme_light() +
  xlab("") + ylab('β-diversity (between pairs of sites)')
  

ggall <- 
  gga / ggb +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
ggall

ggsave(ggall, width = 7, height = 9, dpi = 500,
       filename = paste0("results/2_species_diversity/_S5_violin_plot_alpha_beta_SD_by_site.png"))
ggsave(ggall, width = 7, height = 9, 
       filename = paste0("results/2_species_diversity/_S5_violin_plot_alpha_beta_SD_by_site.pdf"))

