
# read distance matrix ----
level = "site"
comm_delin = "taxonomic_scale_Fishbase"

sd_mat <- readRDS(paste0("results/2_species_diversity/sd_list_pairwise_", level, "_", comm_delin, ".RDS"))



# PCoA ----

gg_list <- list()
i=1
for (metricSD in c("beta.jac", "beta.jtu")){
  # for (comm in names(sd_mat)){
  for (comm in names(sd_mat)[c(2,4)]){
      # read distance matrix
    mat_SDbeta <- sd_mat[[comm]][[metricSD]]
    
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
      scale_color_manual(values = color_perso) +
      ggtitle(comm, subtitle = metricSD) +
      theme(legend.position="none")
    
    i=i+1
  }
  
}

gg_grob <- arrangeGrob(grobs = gg_list, 
                       ncol = length(names(sd_mat))/2, 
                       # ncol = length(names(sd_mat)), 
                       nrow = 2 # nrow = length(names(sd_mat[[1]])),
          )

plot(gg_grob)

ggsave(gg_grob, width = 10, height = 10, 
       filename = paste0("results/2_species_diversity/PCoA_beta_species_diversity_", comm_delin, ".png"))



# NMDS ----

gg_list <- list()
i=1
for (metricSD in c("beta.jac", "beta.jtu")){
  # for (comm in names(sd_mat)){
  for (comm in names(sd_mat)[c(2,4)]){
    # read distance matrix
    mat_SDbeta <- sd_mat[[comm]][[metricSD]]
    
    # run pcoa
    NMDS1 <- ecodist::nmds(mat_SDbeta, mindim = 1, maxdim = 2)
    NMDS2 <- vegan::metaMDS(mat_SDbeta)
    
    # plot
    ecodist::plot.nmds(NMDS1)
    temp <- nmds.min(NMDS1)
    ordiplot (NMDS2)
    stressplot (NMDS2)
    
  }
  
}

