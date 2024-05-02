
# Read data ----
comm_delin <- "taxonomy"
names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")

data_traits <- read_csv("data/traits_Luiz_et_al_2013.csv")

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



# Trait by community ----
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



# Summary ----
pldSUMMARY <- 
  data0_scales %>% 
  group_by(Community) %>% 
  summarise(n = n(),
            mean_PLD = mean(MeanPLD_days),
            sd_PLD = sd(MeanPLD_days),
            mean_Bodysize = mean(Bodysize_cm),
            sd_Bodysize = sd(Bodysize_cm))
pldSUMMARY


# Plot ----
## {FIGURE 3} ####
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





