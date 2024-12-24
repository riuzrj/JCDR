library(rgl)
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
color_vector <- c('normal' = "green", 'ulcerative colitis' = "blue", "Crohn's disease" = "red")

colors_to_use <- color_vector[proj_data$group]

plot3d(proj_data$PC1, proj_data$PC2, proj_data$PC3,
       col=colors_to_use, size=3)

# Projection of PLS
pls <- spp.project.pls(X, Y, r = 3, xfm=FALSE, xfm.opts=list(), robust=FALSE)
proj_data_pls = spp.embed(X, pls$A)
proj_data_pls = cbind(proj_data_pls, group)
proj_data_pls = as.data.frame(proj_data_pls)
# 修改列名
colnames(proj_data_pls) = c('PC1', 'PC2', 'PC3', 'group')
plot3d(proj_data_pls$PC1, proj_data_pls$PC2, proj_data_pls$PC3,
       col=colors_to_use, size=3)


# Projection of PCA+PLS
robust.first <- TRUE
info <- spp.utils.info(X, Y, robust=robust.first)
priors <- info$priors; centroids <- info$centroids
first.moment.proj <-spp.utils.deltas(centroids, priors, robust=robust.first)

# lrlda_1 <- lol.project.lrlda(X, Y, r = 4, xfm=FALSE, xfm.opts=list(), robust=FALSE)
# pls_1 <- lol.project.pls(X, Y, r = 4, xfm=FALSE, xfm.opts=list(), robust=FALSE)
# pca_pls <- cbind(lrlda_1$A[, 2], pls_1$A[, 6])
# pca_pls <- cbind(first.moment.proj, pca_pls)
# proj_data_pca_pls = lol.embed(X, pca_pls)
# proj_data_pca_pls = cbind(proj_data_pca_pls, group)
# proj_data_pca_pls = as.data.frame(proj_data_pca_pls)

lrlda_1 <- spp.project.lrlda(X, Y, r = 30, xfm=FALSE, xfm.opts=list(), robust=FALSE)
proj_data_pca_pls <- cbind(first.moment.proj, lrlda_1$A[, 3])
proj_data_pca_pls = spp.embed(X, proj_data_pca_pls)
proj_data_pca_pls = cbind(proj_data_pca_pls, group)
proj_data_pca_pls = as.data.frame(proj_data_pca_pls)
colnames(proj_data_pca_pls) = c('PC1', 'PC2', 'PC3', 'group')
#
# colors <- c("red", "green", "blue")
# col_vector <- colors[proj_data_pca_pls$group]
# open3d()
# x = proj_data_pca_pls$PC1
# y = proj_data_pca_pls$PC2
# z = proj_data_pca_pls$PC3
# spheres3d(x, y, z, radius=0.05, color=col_vector, alpha=0.5)
#
#

# axes3d(edges = c("bbox", "zneg", "yneg", "xpos"), col="black", labels=c("X", "Y", "Z"))
# title3d(main="3D Scatterplot with Different Colors for Categories", line=2)

color_vector <- c('normal' = "orange", 'ulcerative colitis' ="skyblue", "Crohn's disease" = "purple")

colors_to_use <- color_vector[proj_data_pca_pls$group]
plot3d(proj_data_pca_pls$PC1, proj_data_pca_pls$PC2, proj_data_pca_pls$PC3,
       col=colors_to_use, size = 1.5,type = 's')

