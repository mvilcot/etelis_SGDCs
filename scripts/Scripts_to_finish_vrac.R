# ---- comparative taxonomy !!NOT FINISHED!! ----

data_species2 <- read.csv("data/data_species.csv")
data_fishtree <- read.csv("data/PFC_taxonomy.csv")
data_fishbase <- rfishbase::load_taxa()

Lutjanidae <- list()

Lutjanidae$lutj_agnes <-
  data_species %>%
  dplyr::filter(family == "Lutjanidae") %>%
  pull(species) %>%
  as.factor()

Lutjanidae$lutj_zurich <-
  data_species2 %>%
  dplyr::filter(family == "Lutjanidae") %>%
  pull(species) %>%
  as.factor()

Lutjanidae$lutj_fishtree <-
  data_fishtree %>%
  dplyr::filter(family == "Lutjanidae") %>%
  pull(genus.species) %>%
  sub(pattern = " ", replacement = "_") %>%
  as.factor()

Lutjanidae$lutj_fishbase <-
  data_fishbase %>%
  dplyr::filter(Family == "Lutjanidae") %>%
  pull(Species) %>%
  sub(pattern = " ", replacement = "_") %>%
  as.factor()

Lutj_all <- factor(levels(c(Lutjanidae[[1]], Lutjanidae[[2]], Lutjanidae[[3]], Lutjanidae[[4]])))

temp <- sapply(Lutjanidae, function(x) table(x))
test <- do.call(cbind, temp)
