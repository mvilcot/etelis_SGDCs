# nMDS
#############################################################################################-
# 1. Analyse en Coordonn�es Principales PCoA - metric Multidimensional scaling MDS ####
#############################################################################################-
# Torgerson,  W.  S.  (1952). Multidimensional scaling:  I.  Theory  and  method.  Psychometrika,  17,  401-419. 
# Torgerson,  W.  S.  (1958).  Theory  and  methods of  scaling.  New  York:  Wiley. 
# Gower,  J.  C.  (1966).  Some  distance  properties  of  latent  root  and  vector  methods  used  in multivariate  analysis.  Biometrika,  53,  325-338

# Attention PCoA-MDS est tres diff�rent du nMDS d'un point de vue m�thodologique. 
# Confusion existe parfois dans la litt�rature au niveau des noms ("MDS" alors que c'est un "nMDS" qui a �t� r�alise !)

require(vegan)
require(ade4)

data(dune)
data(dune.env)

# Calcul Bray Curtis (en version metrique pour �tre utilis� dans PCoA). Notes : la plupart des indices non-m�triques utilises en ecologie peuvent le devenir par une transformation racine carr�e (cf CM). Toutefois il est important de verifier a quel point cette transformation peut distordre l'information initiale
# Dissimilarit� de Bray Curtis
bcu=vegdist(dune,"bray")

# Distance de Bray Curtis
bc=sqrt(vegdist(dune,"bray"))

# Comparaison des 2 matrices afin d'evaluer si la transformation racine carree a engendre une distorsion des donnees pour rendre la matrice metrique
plot(bc,bcu,xlim=c(0,1), ylim=c(0,1))
abline(0,1)

# Verifier si propriete euclidienne respectee, car sinon la PCoA peut generer des valeurs propres negatives sur les derniers axes (probleme s'ils doivent etre retenus)
is.euclid(bc)

# si la reponse est False, on applique une transformation de Cailliez (ajout de coefficients aux valeurs initiales)
# bc1=cailliez(bc)

# Puis il aurait fallu v�rifier le degre de distorsion de l'information engendree par la transformation de Cailliez
# plot(bc1,bcu,xlim=c(0,max(bc1,bcu)), ylim=c(0,max(bc1,bcu)))
# abline(0,1)

# plot(bc1,bc,xlim=c(0,max(bc1,bc)), ylim=c(0,max(bc1,bc)))
# abline(0,1)


# PCoA
pcbc=dudi.pco(bc)

# Pourcentage associ� aux axes
pourc=round((pcbc$eig/sum(pcbc$eig))*100,2)
pourc
cumsum(pourc)

# Projections des individus
s.label(pcbc$li, sub="Bray curtis")
# Calcul des distances euclidiennes entre objets sur le 1er plan factoriel
de=dist(pcbc$li[,1:2])

# Verifier que les individus sont bien representes (projetes) � partir de la valeur de leur contribution relative

cont=inertia.dudi(pcbc, row.inertia = TRUE)
cos2=abs(cont$row.rel)/10000
par(mfrow = c(1,2))
plot(pcbc$li[,1], cos2[,1])
plot(pcbc$li[,2], cos2[,2])
# Les objets au centre de chaque axe ont des contributions relatives faibles

# Variables explicatives,  nombre de modalit�s et echantillons par modalite
table(dune.env$Management)
table(dune.env$Use)

# Interaction des 2 facteurs
int=as.factor(paste(dune.env$Management, dune.env$Use, sep="_"))
table(int)  # peu d'observations par modalit�

par(mfrow = c(1,3))
s.class(pcbc$li, dune.env$Management, col=c(1:4))
s.class(pcbc$li, dune.env$Use, col=c(1:3))
s.class(pcbc$li, int, col=c(1:11))

# Anderson, M.J. 2001. A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26: 32�46.  DOI: 10.1111/j.1442-9993.2001.01070.pp.x

deu=dist(pcbc$li[,1:2])
require(vegan)
adonis(deu ~ dune.env$Management, data=dune.env, permutations=999)
adonis(deu ~ dune.env$Use, data=dune.env, permutations=999)
adonis(deu ~ dune.env$Use*dune.env$Management, permutations=999)

# Post Hoc Permanova n'exite pas (c'est la representation sur la PCoA qui permet d'evaluer les differences entre modalites)

# Projection a posteriori de la contribution des variables (correlations des esp�ces aux axes)
require(ape)
pcobc=pcoa(bc)
biplot(pcobc, dune)
# Afin d'�viter de mauvaises interpr�tations du biplot, il serait preferable de faire une repr�sentation separee des objets et variables (pour cela il est possible de reprendre le code source fonction biplot() et de le modifier)

###########-
# 2. nMDS ####
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
# Application 2 ####
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



##############################################################-
# Test de Mantel : corr�lation entre deux matrices de distance ####
##############################################################-



# - Exemple 1 : Relation entre la distance geographique et la distance genetique entre differents villages indiens Yanomama

data(yanomama)
gen <- as.dist(yanomama$gen)
geo <- as.dist(yanomama$geo)
plot(geo,gen)

lm1=lm(gen ~ geo)
abline(lm1$coefficients, col=2)
summary(lm1)

r1<- mantel.randtest(geo,gen, nrepet=999) # uniquement avec le r de pearson 

plot(r1, main = "Mantel's test")

# ou 

mantel(geo,gen, method= � pearson �) avec vegan 

# - Exemple 2 : relation entre ecomorphologie et taxonomie

data(ecomor)
names(ecomor)
distmorpho<-dist(scale(ecomor$morpho))
dtaxo <- dist.taxo(ecomor$taxo)

plot(distmorpho, dtaxo)

mantel(distmorpho, dtaxo, method="pearson", permutations=999)

# - Exemple 3 : relation entre composition des assemblages et les conditions environnementales

data(doubs)

poi<-doubs$fish[-c(8),] 
mil<-doubs$env[-c(8),]

distenv<-dist(scale(mil))
distpoi<-vegdist(poi, method="bray")

plot(distenv, distpoi)

mantel(distpoi, distenv, method="pearson", permutations=999) 

lm3=lm(distpoi ~ distenv)
abline(lm3$coefficients, col=2)
summary(lm3)


# OU SI BESOIN en cas de relation non lineaire et/ou points atypiques
mantel(distpoi, distenv, method="spearman")
