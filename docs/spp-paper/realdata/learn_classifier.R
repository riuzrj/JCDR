library(MASS)
library(e1071)
library(adabag)
library(randomForest)
library(class)
# 尝试运行时不加载 'rgl' 的代码
library(adabag, exclude = "rgl")  # 假设 'adabag' 中的函数不依赖 'rgl'

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

# recombine the data
robust.first <- TRUE
info <- spp.utils.info(X, Y, robust=robust.first)
priors <- info$priors; centroids <- info$centroids
first.moment.proj <- spp.utils.deltas(centroids, priors, robust=robust.first)
PC1 <- spp.embed(X, first.moment.proj[, 1])

data <- data.frame(PC1 = PC1, group = Y)

#划分训练集和测试集
set.seed(123)
n <- nrow(data)
train_indices <- sample(1:n, size = 0.7 * n)  # 70% 数据用于训练
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]
dim(test_data)
# 在训练集上训练 LDA 模型
lda_model <- lda(group ~ PC1, data = train_data)
test_data$predicted_group <- predict(lda_model, newdata = test_data)$class
# 在训练集上训练 svm 模型
svm_model <- svm(group ~ PC1, data = train_data, kernel = "linear")
test_data$predicted_group <- predict(svm_model, newdata = test_data)
# 在训练集上训练 Adaboost 模型
adaboost_model <- adabag::boosting(group ~ PC1, data = train_data, mfinal = 100)
test_data$predicted_group <- predict(adaboost_model, newdata = test_data)$class
# 在训练集上训练 Random Forest 模型
rf_model <- randomForest(group ~ PC1, data = train_data, ntree = 100)
test_data$predicted_group <- predict(rf_model, newdata = test_data)
# 在训练集上训练 KNN 模型
train_data$PC1 <- as.matrix(train_data$PC1)
test_data$PC1 <- as.matrix(test_data$PC1)
knn_model <- knn(train = train_data$PC1, test = test_data$PC1, cl = train_data$group, k = 3)
test_data$predicted_group <- knn_model

ggplot(test_data, aes(x = PC1, fill = predicted_group, color = predicted_group)) +
  geom_density(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("KIRC" = "#003fc8", "LGG" = "#ff6d00")) +
  scale_fill_manual(values = c("KIRC" = "#003fc8", "LGG" = "#ff6d00")) +
  labs(x = "Dimension 1", y = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.line.x = element_line(color = "black"),  # 绘制x轴线
    axis.line.y = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "right"
  )
