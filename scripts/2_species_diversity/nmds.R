
# ---- read distance matrix ----
level = "site"
comm_delin = "taxonomy"

sd_mat <- readRDS(paste0("results/2_species_diversity/sd_list_pairwise_", level, "_", comm_delin, ".RDS"))


# ---- NMDS ----

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





# ---- script récupéré je ne sais plus où ----
###########-
## -- 2. nMDS ----
###########-

require(labdsv)

# nMDS sur matrice initiale de dissimilarite de Bray-Curtis (semi-metrique) calcule dans l'application ci-dessus
bn2=bestnmds(bcu,k=2)
bn4=bestnmds(bcu,k=4)

# Comparaison des resultats fournis par MDS et nMDS, notamment les 2 matrices euclidiennes entre objets qu'elles fournissent sur le 1er plan. 
par(mfrow=c(1,4))

# Plot MDS
plot(pcbc$li[,1:2], type="n", sub="MDS")
text(pcbc$li[,1:2], row.names(dune))
abline(h=0,v=0)

# Plot nMDS
plot(bn2, type="n", sub="nMDS")
text(bn2$points, row.names(dune))

# Calcul des distances euclidiennes entre objets sur le 1er plan 
de1=dist(bn2$points[,1:2])

# Comparaison des resultats entre MDS et nMDS sur le 1er plan, i.e. comparaison des 2 matrices euclidiennes entre objets 
plot(de,de1)
# idem en matrix plot avec 1er bissectrice
plot(de, de1,xlim=c(0,max(de, de1)), ylim=c(0,max(de, de1)))
abline(0,1)

# Comparaison distances euclidiennes MDS puis nMDS sur le 1er plan, par rapport � la matrice initiale de dissimilarit� de Bray Curtis
par(mfrow=c(1,2))

plot(de, bcu,xlim=c(0,max(de, bcu)), ylim=c(0,max(de, bcu)), sub="MDS")
abline(0,1)

plot(de1, bcu,xlim=c(0,max(de1, bcu)), ylim=c(0,max(de1, bcu)), sub="nMDS")
abline(0,1)

# On constate bien que l'objectif du nMDS n'est pas de conserver au mieux les valeurs initiales des distances, mais plutot les rangs inter-objets. Sur ces donnees, on constate que la repr�sentation inter-objets sur le 1er plan est proche entre le MDS et le nMDS (i.e. distance relative entre objets, "groupes" form�s). De fait il est pr�f�rable d'utiliser un MDS qui fournit une solution unique pour une matrice initiale (ici Bray Curtis), et pour lequel les axes correspondent � des axes factoriels auxquels sont associes un % d'inertie (pas le cas du nMDS).
# Il serait interessant de faire la meme analyse sur 4 axes


####################-
## --Application 2 ----
####################-

require(labdsv)
data(bryceveg)

# Calcul de la dissimilarite de Bray Curtis (version semi-metrique).
bc2=vegdist(bryceveg,"bray")

bn2=bestnmds(bc2,k=2)
bn4=bestnmds(bc2,k=4)
# Stress superieur pour k=2 dimensions etant donne que contrainte plus forte � representer le plus fidelement possible les dissimilarites pour k=2 que pour k=4.

par(mfrow=c(2,3))
plot(bn2)
de1=dist(bn2$points[,1:2])
plot(bc2, de1)
plot(bc2, de1,xlim=c(0,max(de1,bc2)), ylim=c(0,max(de1,bc2))) # normal car la conservation des val. initiales des dissimilarites n'est pas l'objecif du nMDS
abline(0,1)

plot(bn4)
de2=dist(bn4$points[,1:4])
plot(bc2, de2)
plot(bc2, de2,xlim=c(0,max(de2,bc2)), ylim=c(0,max(de2,bc2)))
abline(0,1)

# Stress superieur pour k=2 dimensions etant donne que contrainte plus forte � representer le plus fidelement possible les dissimilarites pour k=2 que pour k=4.
# Si l'inverse, relancer bestnmds.
# Le nMDS cherche surtout a conserver le rang des ressemblances, et pas leurs valeurs initiales.

data(brycesite)
par(mfrow=c(1,3))

# nMDS en deux dimensions en illustrant en couleur les quadrats avec une altitude sup�rieure � 8000m. 
plot(bn2)
points(bn2,brycesite$elev>8000)

# nMDS en deux dimensions en illustrant en couleur les quadrats avec une pente inferieure � 10�. 
plot(bn2)
points(bn2,brycesite$slope<10)


# nMDS en deux dimensions en illustrant en couleur les quadrats avec profondeur de sol importante.
plot(bn2)
points(bn2,brycesite$depth=="deep")

# Le facteur principal expliquant la diff�rence de v�g�tation entre les quadrats est l altitude.

