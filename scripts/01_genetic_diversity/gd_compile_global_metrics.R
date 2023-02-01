# This script takes the global metrics calculated in chapter 2 and adds in the
# phylogenetic metrics

# load ----
# genetic metrics
gd_table_global <-
  readRDS("./data/genetic/06_data_global.rds")

# intra-species phylo metrics
phylo_metrics <-
  read_csv("./intermediate/gd_phylo_metrics.csv",
           show_col_types = FALSE)

# merge ----
gd_table_global <-
  left_join(gd_table_global,
            phylo_metrics,
            by = "species") 

# apply exclusions ----
gd_table_global <-
  gd_table_global %>% 
  filter(!species %in% species_exclude)

# export ----
write_csv(gd_table_global,
          "./intermediate/gd_table_global.csv")
