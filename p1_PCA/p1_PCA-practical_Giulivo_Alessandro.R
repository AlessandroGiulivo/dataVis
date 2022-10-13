##########################################
## Practical Session PCA
##########################################
# Follow the instructions in this file. Decide wether to change the code appropriately 
# add new lines or create a new clean file if necessary.

# Fill gaps (...) with appropiate code or answers

# -------------------------------- #
# 1) Load data into variables: 
# -------------------------------- #
## Make use of any or read.delim, read.table, readlines or scan functions.

# First of all, set your working directory if necessary.
setwd("C:/Users/aless/Desktop/Bioinformatics/Erasmus BCN/Data Visualisation/p1_PCA")

# Read data
P1_data = read.delim('P1.tsv.gz', row.names = 1)
P2_data = read.delim('P2.tsv.gz', row.names = 1)

## Read metadata
treat = as.factor(readLines('treatment.txt'))

# -------------------------------- #
# 2) Check data and understand format
# -------------------------------- #
print(head(P1_data)[,1:5])
print(head(P2_data)[,1:5])
print(treat)

print(dim(P1_data))
print(dim(P2_data))
print(table(treat))

# -------------------------------- #
# 3) Adapt the data accordingly and 
#    merge both tables of observations
# -------------------------------- #
# Since we need columns of the tables to represent variables and rows to 
# be samples, individuals, we need to adapt the data of the matrices. 
# Then, we merged experiments into a single data.frame.

# Modify the matrices.
X1 = t(P1_data)
X2 = t(P2_data)

# Merge replicates
X = rbind(X1,X2)

# Check that we transformed the data correctly:
print(dim(X1))
print(dim(X2))
print(dim(X))

# -------------------------------- #
# 4) Get metadata associated
# -------------------------------- #
meta_data <- data.frame(row.names = rownames(X), 
                        experiment=as.factor(rep(c("A","B"), each=96)),
                        treatment=treat)

# -------------------------------- #
# 5) Create PCA: Non-scale => scale=FALSE
# -------------------------------- #
# Create PCA and check plots generated to check results:
pca_res = stats::prcomp(X)

# -------------------------------- #
# 6) Plot PCA: 
# -------------------------------- #
# Plot PCA using several packages (ggfortify, base plot, factoextra) and 
# colour according to treatment and/or experiment.

library(ggfortify)
library(factoextra)

plot(pca_res$x)

# Color according to treatment
plot(pca_res$x, pch=19, col=treat)

# use autoplot
ggplot2::autoplot(pca_res, data = meta_data, colour = "treatment") + theme_classic()

# use factoextra
factoextra::fviz_pca_ind(pca_res, geom = "point", col.ind = meta_data$treatment)

# -------------------------------- #
## Q1. According to the PCA, justify whether you think that there is any 
## batch effect in our samples or not 
# -------------------------------- #
## A1. Yes, the batch effect is evident in the PCA plots, and it is caused by the fact that
## we are analyzing data from two different experiments: the differences in the data are mostly 
## caused by differences in experiment procedure rather than differences in the behaviour of the system

# Color according to experiment
plot(pca_res$x, pch=19, col=meta_data$experiment)
# or
factoextra::fviz_pca_ind(pca_res, geom = "point", col.ind = meta_data$experiment)
# or
ggplot2::autoplot(pca_res, data = meta_data, colour = "experiment") + theme_classic()

# Plot variance explained by each principal component
factoextra::fviz_eig(pca_res, main = "Principal components: Variance explained")
# or
plot(pca_res)

# -------------------------------- #
# 7) Create scaled PCA
# -------------------------------- #
#pca_res_scaled <- stats::prcomp(X, scale=TRUE)

# -------------------------------- #
## Q2. Did it work? Why? 
# -------------------------------- #
## A2. It didn't work because we cannot rescale the data
## when some of the values are equal to zero
# -------------------------------- #
## Remove those genes that have not been expressed 
## at all in our samples and try again.**  

print(dim(X))
X=X[, colSums(X)>0] ## tip: rowSums > 0
print(dim(X))
pca_res_scaled <- stats::prcomp(X, scale=TRUE)

# -------------------------------- #
# 8) Plot scaled PCA by experiment and treatment
# -------------------------------- #
  
## treatment
plot(pca_res_scaled$x, pch=19, col=treat)
ggplot2::autoplot(pca_res_scaled, data = meta_data, colour = "treatment") + theme_classic()

## experiment
plot(pca_res_scaled$x, pch=19, col=meta_data$experiment)
ggplot2::autoplot(pca_res_scaled, data = meta_data, colour = "experiment") + theme_classic()

# -------------------------------- #
## Q3. Has the batch effect observed been corrected?
## If not, is there a clear distinction between batches in any of the axes (PC1 and PC2)?
# -------------------------------- #
## A3. No, the batch effect has not been totally removed,
## there is still a distinction between batches in the
## PC1 axis caused by the differences in the experiment 

# -------------------------------- #
# 9) Plot coefficients
# -------------------------------- #
## We decide to look at the coefficients to understand why.
## Read some info for prcomp and factoextra::get_pca_var and see 
## where is stored this information in prcomp.

?prcomp
?factoextra::get_pca_var

## Plot coefficients results for the first PC.
plot(pca_res_scaled$rotation[,1], pch=19, col=treat)

# non-scale
plot(pca_res$rotation[,1], pch=19, col=treat)

# -------------------------------- #
## Q4. Most coefficients in the y axis are negative. What does it mean for a 
## sample to have a positive value for this PC?
# -------------------------------- #
## A4. Cells having a positive value for this PC are cells with very low or no expression.

# -------------------------------- #
# 10) Explore normalization
# -------------------------------- #
  
# We would check for the total gene expression per cell. 
plot(rowSums(X), type="h", col=meta_data$treatment) ## by treatment
plot(rowSums(X), type="h", col=meta_data$experiment) ## by experiment

# -------------------------------- #
# Q5. Is there any difference between batches?
# Is there any difference inside the same treated/untreated group?
# Do you expect cells with a high value of total expression to have a positive or negative value in PC1?
# -------------------------------- #
## A5. Yes, there is a clear difference in the average expression between batches,
## and there are also differences within the same treated/untreated group
## (the treated cells are in general more expressed than the untreated);
## We expect cells with a high value of total expression to have a positive value in PC1.

# -------------------------------- #
# 11) Check for significance
# -------------------------------- #
  
# Get pvalue for correlation
plot(rowSums(X), pca_res_scaled$x[,1])

library(dplyr)
pca_points <- dplyr::as_tibble(pca_res_scaled$x) %>% bind_cols(meta_data) %>% as.data.frame()
pca_points$exprs <- rowSums(X)
pc1_mod <- lm(PC1 ~ exprs, pca_points)
summary(pc1_mod)

# -------------------------------- #
# Q6. Are PC1 and total gene expression correlated?
# An important number of samples seem to have very low or no expression. 
# What do you think this is due to?
# -------------------------------- #
## A6. Yes, they are correlated, since the P-value is very small.
## The many cells which appear to have very low or no expression
## are probably dead cells.

# -------------------------------- #
# 12) Data normalization
# -------------------------------- #
# In RNA-seq experiments  the data is normalized using different factors such as 
# library size (number of genes) or gene length. We could think of total expression 
# per sample as the factor to use to normalize the data. 

Y = X/rowSums(X)
pca_norm <- stats::prcomp(Y, scale=TRUE)

plot(pca_norm$x, pch=19, col=meta_data$experiment)
plot(pca_norm$x, pch=19, col=meta_data$treatment)

ggplot2::autoplot(pca_norm, data=meta_data, colour="experiment") + theme_classic()
ggplot2::autoplot(pca_norm, data=meta_data, colour="treatment") + theme_classic()


# -------------------------------- #
## Q7. According to the resulting plots, explain what groups can be distinguished 
## and what the main variables represented by the PC1 and PC2 are.
# -------------------------------- #
## A7. Now we have three groups of samples:
## - One on the right with very low or no expression (dead cells);
## - One on the top left, the `untreated` cells;
## - One on the bottom left, the `treated` cells.
## So the main variable represented by the PC1 is the total gene expression (dead vs alive cells),
## while the PC2 represents the differences in expression (treated vs untreated cells).

