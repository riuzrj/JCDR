library(scatterplot3d)
library(datamicroarray)

file_source <- list.files('/Users/ruijuanzhong/JCDR/R')
for (i in 1:length(file_source)) {
  file_source[i] <- paste('/Users/ruijuanzhong/JCDR/R/',file_source[i],sep ="")}
sapply(file_source,source)
file_source <- list.files('/Users/ruijuanzhong/JCDR/docs/lolR')
for (i in 1:length(file_source)) {
  file_source[i] <- paste('/Users/ruijuanzhong/JCDR/docs/lolR/',file_source[i],sep ="")}
sapply(file_source,source)

#data
data('burczynski', package = 'datamicroarray')
X = as.matrix(burczynski$x)
Y = as.factor(burczynski$y)
# Projcetion of PCA
pca <- spp.project.lrlda(X, Y, r = 3, xfm=FALSE, xfm.opts=list(), robust=FALSE)
proj_data = spp.embed(X, pca$A)

group = as.matrix(Y)
proj_data = cbind(proj_data, group)
proj_data = as.data.frame(proj_data)

colnames(proj_data) = c('PC1', 'PC2', 'PC3', 'group')

color_vector <- c('normal' = "orange", 'ulcerative colitis' ="skyblue", "Crohn's disease" = "purple")
colors_to_use <- color_vector[proj_data$group]
scatterplot3d(proj_data$PC1, proj_data$PC2, proj_data$PC3,
              color=colors_to_use, pch=19)


# Projection of PLS
pls <- spp.project.pls(X, Y, r = 3, xfm=FALSE, xfm.opts=list(), robust=FALSE)
proj_data_pls = spp.embed(X, pls$A)
proj_data_pls = cbind(proj_data_pls, group)
proj_data_pls = as.data.frame(proj_data_pls)
colnames(proj_data_pls) = c('PC1', 'PC2', 'PC3', 'group')
scatterplot3d(proj_data_pls$PC1, proj_data_pls$PC2, proj_data_pls$PC3,
              color=colors_to_use, pch=19)


# Projection of PCA+PLS
robust.first <- TRUE
info <- spp.utils.info(X, Y, robust=robust.first)
priors <- info$priors; centroids <- info$centroids
first.moment.proj <- spp.utils.deltas(centroids, priors, robust=robust.first)

lrlda_1 <- spp.project.lrlda(X, Y, r = 4, xfm=FALSE, xfm.opts=list(), robust=FALSE)
pls_1 <- spp.project.pls(X, Y, r = 4, xfm=FALSE, xfm.opts=list(), robust=FALSE)
pca_pls <- cbind(lrlda_1$A[, 2], pls_1$A[, 3])
pca_pls <- cbind(first.moment.proj, pca_pls)
proj_data_pca_pls = spp.embed(X, pca_pls)
proj_data_pca_pls = cbind(proj_data_pca_pls, group)
proj_data_pca_pls = as.data.frame(proj_data_pca_pls)
colnames(proj_data_pca_pls) = c('PC1', 'PC2', 'PC3', 'group')
scatterplot3d(proj_data_pca_pls$PC1, proj_data_pca_pls$PC2, proj_data_pca_pls$PC3,
              color=colors_to_use, pch=19)
