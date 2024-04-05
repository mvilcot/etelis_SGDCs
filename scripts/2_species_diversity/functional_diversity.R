library(tidyverse)
library(ggridges)

# Read data ----
data_traits <- read_csv("data/traits_Luiz_et_al_2013.csv")
data_taxo <- read_csv("data/Taxonomy_Fishbase.csv")
list_communities <- readRDS("intermediate/2_species_diversity/List_community_taxonomy.RDS")

data_traits <- 
  data_traits %>% 
  mutate(Species = paste(Genus, Species, sep = "_")) %>% 
  rename(FAMILY = Family) %>% 
  rename(GENUS = Genus) %>% 
  mutate(FAMILY = str_to_title(FAMILY)) %>% 
  left_join(data_taxo, by = join_by(Species)) %>% 
  relocate(Species, GENUS, Genus, Subfamily, FAMILY, Family, Order, Class, SuperClass) %>% 
  filter(Species %in% list_communities$Teleostei)

data_traits %>% filter(FAMILY != Family)


# Boxplot by family ----
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




# Density by family ----
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


# Histograms ----
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


# Boxplot for this study ----

data1_Class <- 
  data_traits %>% 
  filter(Class == "Teleostei") %>% 
  mutate(Scale = Class)

data2_Order <- 
  data_traits %>% 
  filter(Order == "Eupercaria/misc") %>% 
  mutate(Scale = Order)

data3_Subfam <- 
  data_traits %>% 
  filter(Family == "Lutjanidae") %>% 
  mutate(Scale = Family)

data4_Subfam <- 
  data_traits %>% 
  filter(Subfamily == "Etelinae") %>% 
  mutate(Scale = Subfamily)

data0_scales <- 
  rbind(data1_Class,
        data2_Order,
        data3_Subfam,
        data4_Subfam) %>% 
  dplyr::select(Species, MeanPLD_days, Scale) %>% 
  mutate(Scale = gsub("/misc", "", Scale)) %>% 
  mutate(Scale = factor(Scale, levels = c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")))

# data0_scales %>% 
#   ggplot(aes(MeanPLD_days, Scale)) +
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
  group_by(Scale) %>% 
  summarise(n = n(),
            mean_PLD = mean(MeanPLD_days),
            sd_PLD = sd(MeanPLD_days))

pldSUMMARY %>% write_csv("results/2_species_diversity/Functional_diversity_by_community_PLD.csv")





### {FIGURE 3} ####
temp <- 
  data0_scales %>% 
  left_join(pldSUMMARY, by = "Scale") %>% 
  mutate(Scale = paste0(Scale, "\n", "(n=", n, ")")) 

temp <- 
  temp %>%
  mutate(Scale = factor(Scale, levels = rev(unique(TEST$Scale))))
         
temp %>% 
  ggplot(aes(Scale, MeanPLD_days)) +
  geom_violin(fill = "grey50", color = "grey50") +
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = 45, color = "darkorange", linetype = 2) +
  xlab("") + ylab("PLD (days)") +
  theme_light()

ggsave("results/2_species_diversity/_3_Functional_diversity_by_community_PLD.png", 
       width = 4, height = 4, dpi = 500)
ggsave("results/2_species_diversity/_3_Functional_diversity_by_community_PLD.pdf", 
       width = 4, height = 4)


# TESTS STATS ----
data2_scales <- 
  data0_scales %>% 
  filter(Scale != "Etelinae")

# Homogénéité des variances
bartlett.test(data2_scales$MeanPLD_days, data2_scales$Scale)

library(car)
leveneTest(MeanPLD_days ~ Scale, data = data2_scales)
fligner.test(MeanPLD_days ~ Scale, data = data2_scales)

leveneTest(MeanPLD_days ~ Scale, data = data0_scales %>% filter(Scale %in% c("Teleostei", "Eupercaria")))
leveneTest(MeanPLD_days ~ Scale, data = data0_scales %>% filter(Scale %in% c("Teleostei", "Lutjanidae")))
leveneTest(MeanPLD_days ~ Scale, data = data0_scales %>% filter(Scale %in% c("Eupercaria", "Lutjanidae")))


# ANOVA
model <- lm(data2_scales$MeanPLD_days ~ data2_scales$Scale)
anova(model)

aovc <- aov(MeanPLD_days~Scale, data=data2_scales)
summary(aovc) # le summary() donne ici la table d'anova et non pas les paramètres estimés par le modèle 

# Test de comparaison post-hoc de Tukey
par(mfrow=c(1,1))
hsd <- TukeyHSD(aovc)
hsd
plot(hsd)









