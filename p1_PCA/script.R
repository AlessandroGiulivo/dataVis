P1 <- read.delim('P1.tsv.gz', row.names = 1)
P2 <- read.delim('P2.tsv.gz', row.names = 1)
treat <- as.factor(readLines('treatment.txt'))

dim(P1)
dim(P2)
table(treat)

X1 <- t(P1)
X2 <- t(P2)

X <- rbind(X1, X2)

meta_data <- data.frame(row.names = rownames(X),
                        experiment=as.factor(rep(c("A","B"), each=96)),
                        treatment=treat)

pca = stats::prcomp(X)

plot(pca$x, pch=19, col=treat)

plot(pca$x, pch=19, col=meta_data$experiment)

factoextra::fviz_eig(pca, main = "Principal components: Variance explained")

X=X[, colSums(X)>0]
pca_res_scaled = stats::prcomp(X, scale = TRUE)
print(dim(X))

plot(pca$x, pch=19, col=meta_data$treatment)

pca_res_scaled <- stats::prcomp(X, scale=TRUE)
plot(pca_res_scaled$x, pch=19, col=treat)

plot(pca_res_scaled$x, pch=19, col=meta_data$experiment)
plot(pca_res_scaled$x, pch=19, col=meta_data$treatment)

## plot total gene expression
plot(rowSums(X), type="h", col=treat)
plot(rowSums(X), type="h", col=meta_data$treatment)
plot(rowSums(X), type="h", col=meta_data$experiment)
## Get correlation
plot(rowSums(X), pca_res_scaled$x[,1])
library(dplyr)
library(ggplot2)
pca_points <- dplyr::as_tibble(pca_res_scaled$x) %>% bind_cols(meta_data) %>% as.data.frame()
pca_points$exprs <- rowSums(X)
pc1_mod <- lm(PC1 ~ exprs, pca_points)
summary(pc1_mod)
## Normalize data
Y = X/rowSums(X)
pca_norm <- stats::prcomp(Y, scale=TRUE)
plot(pca_norm$x, pch=19, col=meta_data$experiment)
plot(pca_norm$x, pch=19, col=meta_data$treatment)

ggplot2::autoplot(pca_norm, data=meta_data, colour="experiment") +
  theme_classic()
