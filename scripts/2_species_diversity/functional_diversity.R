
# Read data ----
comm_delin <- "taxonomy"
names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")

data_traits <- read_csv("data/traits_Luiz_et_al_2013.csv")
data_taxo <- read_csv("data/Taxonomy_Fishbase.csv")
# list_communities <- readRDS("intermediate/2_species_diversity/List_community_taxonomy.RDS")
species_present <- 
  readRDS("intermediate/2_species_diversity/PA_Mat_GaspObis_allstations.RDS") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  dplyr::filter(V1 == 1) %>% 
  pull(species)

## join to taxo data
data_traits <- 
  data_traits %>% 
  mutate(Species = paste(Genus, Species, sep = "_")) %>% 
  rename(FAMILY = Family) %>% 
  rename(GENUS = Genus) %>% 
  mutate(FAMILY = str_to_title(FAMILY)) %>% 
  left_join(data_taxo, by = join_by(Species)) %>% 
  filter(Species %in% species_present) %>% #select only species present in the study!!
  relocate(Species, GENUS, Genus, Subfamily, FAMILY, Family, Order, Class, SuperClass)

## check differences beteen data bases
data_traits %>% filter(FAMILY != Family)



# Boxplot for this study ----
data1_Class <- 
  data_traits %>% 
  filter(Class == "Teleostei") %>% 
  mutate(Community = Class)

data2_Order <- 
  data_traits %>% 
  filter(Order == "Eupercaria/misc") %>% 
  mutate(Community = Order)

data3_Family <- 
  data_traits %>% 
  filter(Family == "Lutjanidae") %>% 
  mutate(Community = Family)

data4_Subfam <- 
  data_traits %>% 
  filter(Subfamily == "Etelinae") %>% 
  mutate(Community = Subfamily)

data0_scales <- 
  rbind(data1_Class,
        data2_Order,
        data3_Family,
        data4_Subfam) %>% 
  dplyr::select(Species, MeanPLD_days, Bodysize_cm, Community) %>% 
  mutate(Community = gsub("/misc", "", Community)) %>% 
  mutate(Community = factor(Community, levels = c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")))

# data0_scales %>%
#   ggplot(aes(MeanPLD_days, Community)) +
#   geom_density_ridges(quantile_lines = FALSE, quantiles = 0.5,
#                       jittered_points = TRUE,
#                       position = position_points_jitter(width = 0.05, height = 0),
#                       point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
#   xlab("") +
#   coord_flip() +
#   theme_light()



# Statistics ----
pldSUMMARY <- 
  data0_scales %>% 
  group_by(Community) %>% 
  summarise(n = n(),
            mean_PLD = mean(MeanPLD_days),
            sd_PLD = sd(MeanPLD_days),
            mean_Bodysize = mean(Bodysize_cm),
            sd_Bodysize = sd(Bodysize_cm))
pldSUMMARY

# pldSUMMARY %>% write_csv("results/2_species_diversity/Functional_diversity_by_community_PLD.csv")




## {FIGURE 3} ####
### PLD ----
temp <- 
  data0_scales %>% 
  left_join(pldSUMMARY %>% dplyr::select(Community,n), by = "Community") %>% 
  mutate(Community = paste0(Community, "\n", "(n=", n, ")")) 

temp <- 
  temp %>%
  mutate(Community = factor(Community, levels = rev(unique(temp$Community))))
         
temp %>% 
  ggplot(aes(Community, MeanPLD_days)) +
  geom_violin(fill = "grey50", color = "grey50") +
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = 45, color = "darkorange", linetype = 2) +
  xlab("") + ylab("PLD (days)") +
  theme_light()

ggsave("results/2_species_diversity/_3_Functional_diversity_by_community_PLD_speciespresent.png", 
       width = 4, height = 4, dpi = 500)
ggsave("results/2_species_diversity/_3_Functional_diversity_by_community_PLD_speciespresent.pdf", 
       width = 4, height = 4)


temp %>% 
  pivot_longer(cols = -c(Species, Community, n)) %>%
  ggplot(aes(Community, value)) +
  geom_violin(fill = "grey50", color = "grey50") +
  geom_boxplot(width=0.1) +
  facet_grid2(name ~ ., scales = "free_y") +
  xlab("") + ylab("") +
  theme_light()

ggsave("results/2_species_diversity/_3_Functional_diversity_by_community_PLD_BS_speciespresent.png", 
       width = 8, height = 8, dpi = 500)
ggsave("results/2_species_diversity/_3_Functional_diversity_by_community_PLD_BS_speciespresent.pdf", 
       width = 8, height = 8)






# TESTS STATS ----
data2_scales <- 
  data0_scales %>% 
  filter(Community != "Etelinae")

# Homogénéité des variances
bartlett.test(data2_scales$MeanPLD_days, data2_scales$Community)

library(car)
leveneTest(MeanPLD_days ~ Community, data = data2_scales)
fligner.test(MeanPLD_days ~ Community, data = data2_scales)

leveneTest(MeanPLD_days ~ Community, data = data0_scales %>% filter(Community %in% c("Teleostei", "Eupercaria")))
leveneTest(MeanPLD_days ~ Community, data = data0_scales %>% filter(Community %in% c("Teleostei", "Lutjanidae")))
leveneTest(MeanPLD_days ~ Community, data = data0_scales %>% filter(Community %in% c("Eupercaria", "Lutjanidae")))


# ANOVA
model <- lm(data2_scales$MeanPLD_days ~ data2_scales$Community)
anova(model)

aovc <- aov(MeanPLD_days~Community, data=data2_scales)
summary(aovc) # le summary() donne ici la table d'anova et non pas les paramètres estimés par le modèle 

# Test de comparaison post-hoc de Tukey
par(mfrow=c(1,1))
hsd <- TukeyHSD(aovc)
hsd
# plot(hsd)




# *** Boxplot by family ----
data_traits %>% 
  ggplot(aes(Family, MeanPLD_days)) + 
  geom_boxplot(aes(fill = Family), show.legend = F) +
  facet_wrap(~ Order, scales = "free_x") +
  labs(x="Family", y="PLD (days)") +
  scale_fill_manual(values = c(Lutjanidae = "blue", Haemulidae= "green")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

data_traits %>% 
  filter(Order == "Eupercaria/misc") %>% 
  ggplot(aes(Family, MeanPLD_days)) + 
  geom_boxplot(aes(fill = Family), show.legend = F) +
  facet_wrap(~ Order, scales = "free_x") +
  labs(x="Family", y="PLD (days)") +
  scale_fill_manual(values = c(Lutjanidae = "blue", Haemulidae= "green")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))




# *** Density by family ----
library(ggridges)
data_traits %>% 
  filter(Order == "Eupercaria/misc") %>%
  ggplot(aes(MeanPLD_days, Family)) + 
  geom_density_ridges(aes(fill = Family), show.legend = FALSE,
                      quantile_lines = FALSE, quantiles = 0.5,
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  geom_vline(xintercept = 45) +
  # facet_wrap(~ Order, scales = "free_y") +
  scale_fill_manual(values = c(Lutjanidae = "blue", Haemulidae= "green")) +
  labs(y="Family", x="PLD (days)") +
  theme_light()

data_traits %>% 
  filter(Family == "Lutjanidae") %>% 
  print(n = 30)




# **** Histograms ----
## Subfamily Etelinae ----
gg1 <- 
  data_traits %>% 
  filter(Subfamily == "Etelinae") %>%
  ggplot() + 
  geom_histogram(aes(x = MeanPLD_days), position = "identity") +
  geom_vline(xintercept = 45) +
  xlim(c(0, 85)) +
  labs(x="") +
  theme_light() +
  ggtitle("Etelinae")


## Family Lutjanidae ----
gg2 <- 
  data_traits %>% 
  filter(Family == "Lutjanidae") %>%
  ggplot() + 
  geom_histogram(aes(x = MeanPLD_days), position = "identity") +
  geom_vline(xintercept = 45) +
  xlim(c(0, 85)) +
  labs(x="") +
  theme_light() +
  ggtitle("Lutjanidae")


## Order Eupercaria/misc ----
gg3 <- 
  data_traits %>% 
  filter(Order == "Eupercaria/misc") %>%
  ggplot() + 
  geom_histogram(aes(x = MeanPLD_days), position = "identity") +
  geom_vline(xintercept = 45) +
  xlim(c(0, 85)) +
  labs(x="") +
  theme_light() +
  ggtitle("Eupercaria/misc")

## Class Teleostei ----
gg4 <- 
  data_traits %>% 
  filter(Class == "Teleostei") %>%
  ggplot() + 
  geom_histogram(aes(x = MeanPLD_days), position = "identity") +
  geom_vline(xintercept = 45) +
  xlim(c(0, 85)) +
  labs(x="PLD (days)") +
  theme_light() +
  ggtitle("Teleostei")


gg1 / gg2 / gg3 / gg4 






