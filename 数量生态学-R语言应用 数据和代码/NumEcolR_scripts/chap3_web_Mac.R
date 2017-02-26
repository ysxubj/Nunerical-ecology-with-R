### CHAPTER 3: ASSOCIATION MEASURES ###

# Load required libraries
library(ade4)
library(vegan)	# should be loaded after ade4 to avoid some conflicts
library(gclus)
library(cluster)
library(FD)

# Import the data from CSV files
spe <- read.csv("DoubsSpe.csv", row.names=1)
env <- read.csv("DoubsEnv.csv", row.names=1)
spa <- read.csv("DoubsSpa.csv", row.names=1)
# Remove empty site 8
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]


# Dissimilarity and distance measures for (semi-)quantitative data
# ****************************************************************

# Bray-Curtis dissimilarity matrix on raw species data
spe.db <- vegdist(spe)	# Bray-Curtis dissimilarity (default)
head(spe.db)
# Bray-Curtis dissimilarity matrix on log-transformed abundances
spe.dbln <- vegdist(log1p(spe))
head(spe.dbln)
# Chord distance matrix
spe.norm <- decostand(spe, "nor")
spe.dc <- dist(spe.norm)
head(spe.dc)
# Hellinger distance matrix
spe.hel <- decostand(spe, "hel")
spe.dh <- dist(spe.hel)
head(spe.dh)


# Dissimilarity measures for binary data
# **************************************

# Jaccard dissimilarity matrix using function vegdist()
spe.dj <- vegdist(spe, "jac", binary=TRUE)
head(spe.dj)
head(sqrt(spe.dj))
# Jaccard dissimilarity matrix using function dist()
spe.dj2 <- dist(spe, "binary")
head(spe.dj2)
# Jaccard dissimilarity matrix using function dist.binary()
spe.dj3 <- dist.binary(spe, method=1)
head(spe.dj3)
# Sorensen dissimilarity matrix using function dist.binary()
spe.ds <- dist.binary(spe, method=5)
head(spe.ds)
# Sorensen dissimilarity matrix using function vegdist()
spe.ds2 <- vegdist(spe, binary=TRUE)
head(spe.ds2)
head(sqrt(spe.ds2))
# Ochiai dissimilarity matrix
spe.och <- dist.binary(spe, method=7)
head(spe.och)


# Graphical display of association matrices
# *****************************************

# The gclus package is required and may be called now, although
# it is called internally by coldiss()
library(gclus)

# Colour plots (also called heat maps, or trellis diagrams in the data 
# analysis literature) using the coldiss() function
# ********************************************************************

# Source the coldiss() function 
source("coldiss.R")  # If necessary, add the path to the file

# Bray-Curtis dissimilarity matrix (on raw data)
# 4 colours with equal-length intervals (useful for comparisons)
quartz(title="Bray-Curtis (raw data)",10,5)
coldiss(spe.db, byrank=FALSE, diag=TRUE)

# Same but on log-transformed data
quartz(title="Bray-Curtis [ln(y+1) data]",10,5)
coldiss(spe.dbln, byrank=FALSE, diag=TRUE)

# Chord distance matrix
quartz(title="Chord",10,5)
coldiss(spe.dc, byrank=FALSE, diag=TRUE)

# Hellinger distance matrix
quartz(title="Hellinger",10,5)
coldiss(spe.dh, byrank=FALSE, diag=TRUE)

# Jaccard distance matrix
quartz(title="Jaccard",10,5)
coldiss(spe.dj, byrank=FALSE, diag=TRUE)

# Simple matching dissimilarity
# (called the Sokal and Michener index in ade4)
spe.s1 <- dist.binary(spe, method=2)
quartz(title="S1 on species data",10,5) 
coldiss(spe.s1^2, byrank=FALSE, diag=TRUE)



# Remove the 'das' variable from the env dataset
env2 <- env[,-1]

# Euclidean distance matrix of the standardized env2 data frame
env.de <- dist(scale(env2))
quartz(title="Environment",10,5)
coldiss(env.de, diag=TRUE)

# Hellinger distance matrix of the species data (equal-sized categories)
quartz(title="Species",10,5)
coldiss(spe.dh, diag=TRUE)

# Euclidean distance matrix on spatial coordinates (2D)
spa.de <- dist(spa)
quartz(title="x-y",10,5)
coldiss(spa.de, diag=TRUE)

# Euclidean distance matrix on distance from the source (1D)
das.df <- as.data.frame(env$das, row.names=rownames(env))
riv.de <- dist(das.df)
quartz(title="Distance from source",10,5) 
coldiss(riv.de, diag=TRUE)



# Compute five binary variables with 30 objects each. Each variable
# has a predefined number of 0 and 1
# Variable 1: 10 x 1 and 20 x 0; the order is randomized
var1 <- sample(c(rep(1,10), rep(0,20)))
# Variable 2: 15 x 0 and 15 x 1, one block each
var2 <- c(rep(0,15), rep(1,15))
# Variable 3: alternation of 3 x 1 and 3 x 0 up to 30 objects
var3 <- rep(c(1,1,1,0,0,0),5)
# Variable 4: alternation of 5 x 1 and 10 x 0 up to 30 objects
var4 <- rep(c(rep(1,5), rep(0,10)), 2)
# Variable 5: 16 objects with randomized distribution of 7 x 1 
# and 9 x 0, followed by 4 x 0 and 10 x 1
var5.1 <- sample(c(rep(1,7), rep(0,9)))
var5.2 <- c(rep(0,4), rep(1,10))
var5 <- c(var5.1, var5.2)

# Variables 1 to 5 are put into a data frame
dat <- data.frame(var1, var2, var3, var4, var5)
dim(dat)

# Computation of a matrix of simple matching coefficients
# (called Sokal and Michener index in ade4)
dat.s1 <- dist.binary(dat, method=2)
quartz(title="S1 on fictitious data",10,5) 
coldiss(dat.s1, diag=TRUE)


# Fictitious data for Gower (S15) index
# Random normal deviates with zero mean and unit standard deviation
var.g1 <- rnorm(30,0,1)
# Random uniform deviates from 0 to 5
var.g2 <- runif(30,0,5)
# Factor with 3 levels (10 objects each)
var.g3 <- gl(3,10)
# Factor with 2 levels, orthogonal to var.g3
var.g4 <- gl(2,5,30)
dat2 <- data.frame(var.g1,var.g2,var.g3,var.g4)
summary(dat2)

# Computation of a matrix of Gower dissimilarity using function daisy()

# Complete data matrix (4 variables)
dat2.S15 <- daisy(dat2, "gower")
range(dat2.S15)
quartz(title="S15 on fictitious data - daisy",10,5) 
coldiss(dat2.S15, diag=TRUE)

# Data matrix with the two orthogonal factors only
dat2partial.S15 <- daisy(dat2[,3:4], "gower")
quartz(title="S15 on fictitious data, 2 factors - daisy",10,5) 
coldiss(dat2partial.S15, diag=TRUE)

# What are the dissimilarity values in the dat2partial.S15 matrix?
levels(factor(dat2partial.S15))

# Computation of a matrix of Gower dissimilarity using function gowdis() 
# of package FD
library(FD) # If not already loaded
?gowdis
dat2.S15.2 <- gowdis(dat2)
range(dat2.S15.2)
quartz(title="S15 on fictitious data - gowdis",10,5) 
coldiss(dat2.S15.2, diag=TRUE)

# Data matrix with the two orthogonal factors only
dat2partial.S15.2 <- gowdis(dat2[,3:4])
quartz(title="S15 on fictitious data, 2 factors - gowdis",10,5) 
coldiss(dat2partial.S15.2, diag=TRUE)

# What are the dissimilarity values in the dat2partial.S15.2 matrix? 
levels(factor(dat2partial.S15.2))



# R-mode dissimilarity matrix
# ***************************

# Transpose matrix of species abundances
spe.t <- t(spe)

# Chi-square pre-transformation followed by Euclidean distance
spe.t.chi <- decostand(spe.t, "chi.square")
spe.t.D16 <- dist(spe.t.chi)
quartz(title="D16 on fish species (R-mode)",10,5)
coldiss(spe.t.D16, diag=TRUE)
  
# Jaccard index on fish presence-absence
spe.t.S7 <- vegdist(spe.t, "jaccard", binary=TRUE)
quartz(title="S7 on fish species (R-mode)",10,5) 
coldiss(spe.t.S7, diag=TRUE)

# Pearson r linear correlation among environmental variables
env.pearson <- cor(env)	# default method = "pearson"
round(env.pearson, 2)

# Reorder the variables prior to plotting
env.o <- order.single(env.pearson)

# pairs: a function to plot a matrix of bivariate scatter diagrams
# and correlation coefficients.
# Correlations are given in the upper panel (with significance levels)
source("panelutils.R") # If necessary give path
quartz(title="Linear correlation matrix",10,10)
op <- par(mfrow=c(1,1), pty="s")
pairs(env[,env.o], lower.panel=panel.smooth, upper.panel=panel.cor,
	diag.panel=panel.hist, main="Pearson Correlation Matrix")
par(op)

# Kendall tau rank correlation among environmental variables
env.ken <- cor(env, method="kendall")
env.o <- order.single(env.ken)
quartz(title="Rank correlation matrix",10,10)
op <- par(mfrow=c(1,1), pty="s")
pairs(env[,env.o], lower.panel=panel.smooth, upper.panel=panel.cor,
	method="kendall", diag.panel=panel.hist, main="Kendall Correlation Matrix")
par(op)
