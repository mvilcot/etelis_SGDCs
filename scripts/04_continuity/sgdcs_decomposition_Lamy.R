
# Appendix S2
#
# R function to compute a decomposition of species-genetic diversity correlations (SGDCs) into the contribution of different variables

SGDC.decomp <- function(SD, GD, FACTOR, plot=TRUE, verbose=TRUE)
#
# Description --
# 
# This method decomposes the covariance between species diversity (SD) and genetic diversity (GD) to identify 
# the contributions of different habitat characteristics to the SGDC.
# These contributions include direct contributions of site characteristics as well as the contributions of 
# covariance between pairs of habitat characteristics
# 
# if there is only one environmental variable you can add it as a numeric vector.
# In such a case, the method compute the contribution of the single variables as the product of 
# the standardized regression coefficients of that variable on both SD and GD. Note that in single linear regression,
# standardized regression coefficients equal correlation coefficients between the response and the explanatory variable.
# The case with only one variable is different from a partial correlation.
#
# Recommendation --
#
# This method relies on linear regressions. It is to the user to check that his data
# fits the assumptions of these models (distribution of response variables SD and GD and of model residuals)
#
# Arguments --
# 
# SD: an vector of species diversity (SD)
# GD: an vector of genetic diversity (GD)
# FACTOR: a data.frame containing the environmental variables, or a vector containing only one environmental variable
# This data.frame should not contain any additional variables, as the SGDC decomposition 
# will be parsed to all variables in FACTOR.
# plot=TRUE, will plot the summary 
#
# Author: Thomas Lamy, September 2014
#
# References --
# Lamy et al (2013) Variation in habitat connectivity generates positive correlations between species and genetic diversity in a metacommunity. Molecular ecology, 22: 4445–4456
# Lamy, Laroche et al. in prep. The contribution of species-genetic diversity correlations to the understanding of community assembly rules.
{
require(vegan)
if(class(SD)!="numeric") stop("SD has to be a numeric vector")
if(class(GD)!="numeric") stop("GD has to be a numeric vector")
FACTOR <- as.data.frame(FACTOR)
if(sum(is.na(SD))!=0|sum(is.na(GD))!=0) stop("Missing data are present and should be removed")
if(length(SD)!=length(GD)) stop("Numbers of rows in SD and GD are different")
if(length(SD)!=dim(FACTOR)[1]) stop("Numbers of rows in SD and GD does not match the number of rows in FACTOR")

# compute standardized (i.e. zero mean and unity sd) variables
SD <- decostand(SD, method="standardize")[,1]
GD <- decostand(GD, method="standardize")[,1]
FACTOR <- decostand(FACTOR, method="standardize")

# CASE 1: simple linear regressions with only one environmental variable
if(dim(FACTOR)[2]==1){
	slr.SD <- lm(SD~FACTOR[,1])
	slr.GD <- lm(GD~FACTOR[,1])
	shapi.SD <- shapiro.test(residuals(slr.SD))
	shapi.GD <- shapiro.test(residuals(slr.GD))
	if(verbose) ({
		cat("Shapiro-Wilk normality test on the linear model residuals",'\n')
		cat("Linear model on SD, p-value =",round(shapi.SD$p.value,5),'\n')
		cat("Linear model on GD, p-value =",round(shapi.GD$p.value,5),'\n')
		})
	# compute the SGDC 
	SGDC <- cor(SD, GD)
	# the contribution of the single variable to the SGDC is the product of the 
	# standardized regression coefficients of that variable on both SD and GD
	contri.FACTOR <- slr.SD$coefficients[-1] * slr.GD$coefficients[-1]
	# unexplained part of the SGDC
	SGDC.res <- SGDC - contri.FACTOR
	# is this unexplained part significant?
	toto <- cor.test(residuals(slr.SD), residuals(slr.GD)) 
	if(verbose) cat("Residual SGDC p-value is", round(toto$p.value,5))
	# summary table
	Contribution <- c(contri.FACTOR, SGDC.res, SGDC)
	Contribution <- as.data.frame(Contribution)
	rownames(Contribution) <- c("Var(X)", "Residuals", "SGDC")
	Contribution$Percent <- round(Contribution$Contribution/SGDC*100, 2)

	if(plot){
		## Graphic 
		barplot(Contribution$Percent, horiz = TRUE, col=c("red", "grey", "black"),
		main = "SGDC decomposition", xlab = "% of the SGDC", xlim = c(min(Contribution$Percent)-10,max(Contribution$Percent)+10), cex.names = rownames(Contribution))
	}
}

# CASE 2: multiple linear regressions with two or more environmental variables
# multiple linear regression
if(dim(FACTOR)[2]>=2){
	mlr.SD <- lm(SD~., data=FACTOR)
	mlr.GD <- lm(GD~., data=FACTOR)
	shapi.SD <- shapiro.test(residuals(mlr.SD))
	shapi.GD <- shapiro.test(residuals(mlr.GD))
	if(verbose) ({
		cat("Shapiro-Wilk normality test on the linear model residuals",'\n')
		cat("Linear model on SD, p-value =",round(shapi.SD$p.value,5),'\n')
		cat("Linear model on GD, p-value =",round(shapi.GD$p.value,5),'\n')
	})
	# standardized regression coefficients of variables on SD
	coef.SD <- mlr.SD$coefficients[-1]
	# standardized regression coefficients of variables on GD
	coef.GD <- mlr.GD$coefficients[-1]
	# names of the variables
	var.names <- names(coef.SD)
	# number of variables
	var.nb <- length(names(coef.SD))
	# compute the variance-covariance matrix among environmental variables.
	# as environmental variables have been standardized, this matrix is the 
	# correlation matrix
	mat.cov <- cov(FACTOR)
	# Variance terms
	coef.var <- as.numeric(coef.SD * coef.GD)
	# Covariance terms
	coef.cov1 <- coef.SD %*% t(coef.GD)
	coef.cov2 <- coef.GD %*% t(coef.SD)
	coef.cov <- (coef.cov1 + coef.cov2) * mat.cov[var.names,var.names]
	# formating results
	names.cov = c()
	for(i in 1:(var.nb-1)){
	  for(j in (i+1):var.nb){
    	names.cov <- c(names.cov,   
           			paste("Cov(",var.names[i],", ",var.names[j],")", sep=""))
  		}
	}
	names.var <- paste("Var(",var.names,")", sep="")

	# compute the SGDC 
	SGDC <- cor(SD, GD)
	# unexplained part of the SGDC
	SGDC.res <- SGDC - sum(coef.var) - sum(as.vector(as.dist(coef.cov)))
	# is this unexplained part significant?
	toto <- cor.test(residuals(mlr.SD), residuals(mlr.GD)) 
	if(verbose) cat("Residual SGDC p-value is", round(toto$p.value,5))

	# summary table
	Contribution <- c(coef.var, as.vector(as.dist(coef.cov)), SGDC.res, SGDC)
	Contribution <- as.data.frame(Contribution)
	rownames(Contribution) <- c(names.var, names.cov, "Residuals", "SGDC")
	Contribution$Percent <- round(Contribution$Contribution/SGDC*100, 2)

	if(plot){
		## Graphic 
	  barplot(Contribution$Percent, horiz = TRUE, col=c(rep("red", var.nb), rep("orange", 6), "grey", "black"),
	          main = "SGDC decomposition", xlab = "% of the SGDC", xlim = c(min(Contribution$Percent)-10,max(Contribution$Percent)+10), 
	          names.arg = rownames(Contribution), cex.names=1, las=1)
	}
}
res <- list(Contribution=Contribution, SGDC.res=round(toto$p.value,5))
res
}

