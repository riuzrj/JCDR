library(ggplot2)
library(reshape2)
library(pheatmap)


file_source <- list.files('/Users/ruijuanzhong/JCDR/R')
for (i in 1:length(file_source)) {
  file_source[i] <- paste('/Users/ruijuanzhong/JCDR/R/',file_source[i],sep ="")}
sapply(file_source,source)

query <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_9 <- GDCprepare(query = query, save = F)
dataPrep_9 <- TCGAanalyze_Preprocessing(object = dataPrep1_9,
                                        cor.cut = 0.6)
dataNorm_9 <- TCGAanalyze_Normalization(tabDF = dataPrep_9,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_9 <- TCGAanalyze_Filtering(tabDF = dataNorm_9,
                                    method = "quantile",
                                    qnt.cut =  0.25)
sampleTypes_9 <- c(rep('KIRC', ncol(dataFilt_9)))

# require data
query <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_10 <- GDCprepare(query = query, save = F)
dataPrep_10 <- TCGAanalyze_Preprocessing(object = dataPrep1_10,
                                         cor.cut = 0.6)
dataNorm_10 <- TCGAanalyze_Normalization(tabDF = dataPrep_10,
                                         geneInfo = geneInfoHT,
                                         method = "geneLength")
dataFilt_10 <- TCGAanalyze_Filtering(tabDF = dataNorm_10,
                                     method = "quantile",
                                     qnt.cut =  0.25)
dim(dataFilt_10)
sampleTypes_10 <- c(rep('LGG', ncol(dataFilt_10)))

common_cols_2 <- Reduce(intersect, list(colnames(t(dataFilt_9)), colnames(t(dataFilt_10))))
KIRC <- t(dataFilt_9)[, common_cols_2]
LGG <- t(dataFilt_10)[, common_cols_2]

X <- rbind(KIRC, LGG)
Y <- factor(c(sampleTypes_9, sampleTypes_10))
dim(X)
table(Y)

# learn projection
X_log_transformed <- log2(X + 1)  # 加1避免log(0)
X <- scale(X_log_transformed)
head(X)
robust.first <- TRUE
info <- spp.utils.info(X, Y, robust=robust.first)
priors <- info$priors; centroids <- info$centroids
first.moment.proj <- spp.utils.deltas(centroids, priors, robust=robust.first)
lrlda <- spp.project.lrlda(X, Y, r = 4, xfm=FALSE, xfm.opts=list(), robust=FALSE)
pls <- spp.project.pls(X, Y, r = 4, xfm=FALSE, xfm.opts=list(), robust=FALSE)
pca_pls_loadings <- cbind(pls$loadings[, 1:2], lrlda$A[, 1:2])
colnames(pca_pls_loadings) <- c("Comp 1", "Comp 2", "Comp 3", "Comp 4")
loadings_melted <- melt(pca_pls_loadings)
summary(loadings_melted$value)
print(colnames(loadings_melted))
colnames(loadings_melted) <- c("gene", "component", "value")
# ggplot(loadings_melted, aes(x = gene, y = component, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                        midpoint = 0, limit = c(-0.02, 0.015)) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal()+
#   theme(axis.text.x = element_blank()) +
#   facet_grid(component ~ ., scales = "free_y", space = "free_y")


ggplot(loadings_melted, aes(x = gene, y = component, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.02, 0.015),
                       name = "Loadings Value") +  # 为颜色条添加标题
  labs(x = NULL, y = NULL) +  # 添加图表标题
  theme_minimal(base_size = 12) +  # 使用较小的基准字体
  theme(axis.text.x = element_blank(),  # 移除 X 轴标签
        axis.text.y = element_text(face = "bold"),  # 使 Y 轴标签更显眼
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 图表标题居中
        legend.title = element_text(size = 10),  # 调整颜色条标题大小
        legend.position = "right",
        strip.text.y = element_blank()) +  # 确保颜色条在图表右侧
  facet_grid(component ~ ., scales = "free_y", space = "free_y")  # 仅在 y 轴分面
