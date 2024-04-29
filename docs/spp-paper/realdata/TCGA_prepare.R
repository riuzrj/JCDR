library(SummarizedExperiment)
library(TCGAbiolinks)
library(EDASeq)
library(sesameData)
library(minfi)
library(limma)
library(MultiAssayExperiment)

query <- GDCquery(
  project = "TCGA-CHOL",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1 <- GDCprepare(query = query, save = F)
#sampleTypes_1 <- colData(dataPrep1)$sample_type
dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1,
                                      cor.cut = 0.6)
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "geneLength")
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.25)
sampleTypes <- c(rep('CHOL', ncol(dataFilt)))

query <- GDCquery(
  project = "TCGA-DLBC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_1 <- GDCprepare(query = query, save = F)
#sampleTypes_1 <- colData(dataPrep1)$sample_type
dataPrep_1 <- TCGAanalyze_Preprocessing(object = dataPrep1_1,
                                        cor.cut = 0.6)
dataNorm_1 <- TCGAanalyze_Normalization(tabDF = dataPrep_1,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_1 <- TCGAanalyze_Filtering(tabDF = dataNorm_1,
                                    method = "quantile",
                                    qnt.cut =  0.25)
sampleTypes_1 <- c(rep('DLBC', ncol(dataFilt_1)))

query <- GDCquery(
  project = "TCGA-KICH",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_2 <- GDCprepare(query = query, save = F)
dataPrep_2 <- TCGAanalyze_Preprocessing(object = dataPrep1_2,
                                        cor.cut = 0.6)
dataNorm_2 <- TCGAanalyze_Normalization(tabDF = dataPrep_2,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_2 <- TCGAanalyze_Filtering(tabDF = dataNorm_2,
                                    method = "quantile",
                                    qnt.cut =  0.25)
sampleTypes_2 <- c(rep('KICH', ncol(dataFilt_2)))

query <- GDCquery(
  project = "TCGA-MESO",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_3 <- GDCprepare(query = query, save = F)
dataPrep_3 <- TCGAanalyze_Preprocessing(object = dataPrep1_3,
                                        cor.cut = 0.6)
dataNorm_3 <- TCGAanalyze_Normalization(tabDF = dataPrep_3,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_3 <- TCGAanalyze_Filtering(tabDF = dataNorm_3,
                                    method = "quantile",
                                    qnt.cut =  0.25)
sampleTypes_3 <- c(rep('MESO', ncol(dataFilt_3)))

query <- GDCquery(
  project = "TCGA-READ",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_4 <- GDCprepare(query = query, save = F)
dataPrep_4 <- TCGAanalyze_Preprocessing(object = dataPrep1_4,
                                        cor.cut = 0.6)
dataNorm_4 <- TCGAanalyze_Normalization(tabDF = dataPrep_4,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_4 <- TCGAanalyze_Filtering(tabDF = dataNorm_4,
                                    method = "quantile",
                                    qnt.cut =  0.25)
sampleTypes_4 <- c(rep('READ', ncol(dataFilt_4)))

query <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_5 <- GDCprepare(query = query, save = F)
dataPrep_5 <- TCGAanalyze_Preprocessing(object = dataPrep1_5,
                                        cor.cut = 0.6)
dataNorm_5 <- TCGAanalyze_Normalization(tabDF = dataPrep_5,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_5 <- TCGAanalyze_Filtering(tabDF = dataNorm_5,
                                    method = "quantile",
                                    qnt.cut =  0.25)
sampleTypes_5 <- c(rep('LIHC', ncol(dataFilt_5)))

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_6 <- GDCprepare(query = query, save = F)
dataPrep_6 <- TCGAanalyze_Preprocessing(object = dataPrep1_6,
                                        cor.cut = 0.6)
dataNorm_6 <- TCGAanalyze_Normalization(tabDF = dataPrep_6,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_6 <- TCGAanalyze_Filtering(tabDF = dataNorm_6,
                                    method = "quantile",
                                    qnt.cut =  0.25)
sampleTypes_6 <- c(rep('LUAD', ncol(dataFilt_6)))

query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_7 <- GDCprepare(query = query, save = F)
dataPrep_7 <- TCGAanalyze_Preprocessing(object = dataPrep1_7,
                                        cor.cut = 0.6)
dataNorm_7 <- TCGAanalyze_Normalization(tabDF = dataPrep_7,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_7 <- TCGAanalyze_Filtering(tabDF = dataNorm_7,
                                    method = "quantile",
                                    qnt.cut =  0.25)
sampleTypes_7 <- c(rep('PAAD', ncol(dataFilt_7)))

query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_8 <- GDCprepare(query = query, save = F)
dataPrep_8 <- TCGAanalyze_Preprocessing(object = dataPrep1_8,
                                        cor.cut = 0.6)
dataNorm_8 <- TCGAanalyze_Normalization(tabDF = dataPrep_8,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_8 <- TCGAanalyze_Filtering(tabDF = dataNorm_8,
                                    method = "quantile",
                                    qnt.cut =  0.25)
sampleTypes_8 <- c(rep('PAAD', ncol(dataFilt_8)))

query <- GDCquery(
  project = "TCGA-READ",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value")
GDCdownload(query = query)
dataPrep1_9 <- GDCprepare(query = query, save = F)
betaMatrix <- assays(dataPrep1_9)[[1]]
betaMatrixClean <- na.omit(betaMatrix)
cpgs <- rowSds(as.matrix(betaMatrixClean))
quantileCutoff <- quantile(cpgs, probs = 0.75)
informativeCpGs <- names(cpgs)[cpgs > quantileCutoff]
featureMatrix <- betaMatrixClean[informativeCpGs, ]
dim(featureMatrix)
normFeatureMatrix <- normalizeBetweenArrays(as.matrix(featureMatrix), method = "quantile")
DNA_READ = t(normFeatureMatrix)


query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value")
GDCdownload(query = query)
dataPrep1_10 <- GDCprepare(query = query, save = F)
betaMatrix_1 <- assays(dataPrep1_10)[[1]]
betaMatrixClean_1 <- na.omit(betaMatrix_1)
cpgs_1 <- rowSds(as.matrix(betaMatrixClean_1))
quantileCutoff_1 <- quantile(cpgs_1, probs = 0.75)
informativeCpGs_1 <- names(cpgs_1)[cpgs_1 > quantileCutoff_1]
featureMatrix_1 <- betaMatrixClean_1[informativeCpGs_1, ]
dim(featureMatrix_1)
normFeatureMatrix_1 <- normalizeBetweenArrays(as.matrix(featureMatrix_1), method = "quantile")
DNA_COAD = t(normFeatureMatrix_1 )

CHOL <- t(dataFilt)
DLBC <- t(dataFilt_1)
KICH <- t(dataFilt_2)
MESO <- t(dataFilt_3)
READ <- t(dataFilt_4)
LIHC <- t(dataFilt_5)
LUAD <- t(dataFilt_6)
PAAD <- t(dataFilt_7)
COAD <- t(dataFilt_8)
common_cols <- Reduce(intersect, list(colnames(COAD), colnames(READ)))
CHOL <- CHOL[, common_cols]
DLBC <- DLBC[, common_cols]
KICH <- KICH[, common_cols]
MESO <- MESO[, common_cols]
READ <- READ[, common_cols]
LIHC <- LIHC[, common_cols]
LUAD <- LUAD[, common_cols]
PAAD <- PAAD[, common_cols]
COAD <- COAD[, common_cols]
X <- rbind(COAD, READ)
dim(X)
Y <- factor(c(sampleTypes_4, sampleTypes_8))
