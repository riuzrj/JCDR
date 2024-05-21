library(GGally)
library(datamicroarray)

file_source <- list.files('/Users/ruijuanzhong/spp')
for (i in 1:length(file_source)) {
  file_source[i] <- paste('/Users/ruijuanzhong/spp/',file_source[i],sep ="")}
sapply(file_source,source)

#data
data('burczynski', package = 'datamicroarray')
X = as.matrix(burczynski$x)
Y = as.factor(burczynski$y)
# Projcetion of PCA
pca <- lol.project.lrlda(X, Y, r = 3, xfm=FALSE, xfm.opts=list(), robust=FALSE)
proj_data = lol.embed(X, pca$A)
proj_data = as.data.frame(proj_data)
colnames(proj_data) = c('PC1', 'PC2', 'PC3')
# Make sure your data is in a data frame with PC1, PC2, and PC3 as columns
df <- data.frame(PC1 = proj_data$PC1, PC2 = proj_data$PC2, PC3 = proj_data$PC3)

# Create the pairs plot
pairs(df, main = "Scatter plot matrix with the first three principal components",
      col = "blue")


#Projection of PLS
pls <- lol.project.pls(X, Y, r = 3, xfm=FALSE, xfm.opts=list(), robust=FALSE)
proj_data_pls = lol.embed(X, pls$A)

proj_data_pls = as.data.frame(proj_data_pls)
# 修改列名
colnames(proj_data_pls) = c('PC1', 'PC2', 'PC3')
# Make sure your data is in a data frame with PC1, PC2, and PC3 as columns
df_1 <- data.frame(PC1 = proj_data_pls$PC1, PC2 = proj_data_pls$PC2, PC3 = proj_data_pls$PC3)
# Create the pairs plot
pairs(df_1, main = "Scatter plot matrix with the first three principal components",
      col = "blue")


# Projection of PCA+PLS
robust.first <- TRUE
info <- lol.utils.info(X, Y, robust=robust.first)
priors <- info$priors; centroids <- info$centroids
first.moment.proj <- lol.utils.deltas(centroids, priors, robust=robust.first)

lrlda_1 <- lol.project.lrlda(X, Y, r = 30, xfm=FALSE, xfm.opts=list(), robust=FALSE)
pls_1 <- lol.project.pls(X, Y, r = 30, xfm=FALSE, xfm.opts=list(), robust=FALSE)
pca_pls <- cbind(lrlda_1$A[, 1:2], pls_1$A[, 3:4])
pca_pls <- cbind(first.moment.proj, pca_pls)
proj_data_pca_pls = lol.embed(X, pca_pls)

proj_data_pca_pls = as.data.frame(proj_data_pca_pls)
# 修改列名
colnames(proj_data_pca_pls) = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6')
# Make sure your data is in a data frame with PC1, PC2, and PC3 as columns
df_2 <- data.frame(PC1 = proj_data_pca_pls$PC1, PC2 = proj_data_pca_pls$PC2, PC3 = proj_data_pca_pls$PC3,
                   PC4 = proj_data_pca_pls$PC4, PC5 = proj_data_pca_pls$PC5, PC6 = proj_data_pca_pls$PC6)
# Create the pairs plot
pairs(df_2, main = "Scatter plot matrix with the first three principal components",
      col = "red")
