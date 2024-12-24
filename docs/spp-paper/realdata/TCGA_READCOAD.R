library(SummarizedExperiment)
library(TCGAbiolinks)
library(EDASeq)
library(sesameData)
library(minfi)
library(limma)
library(MultiAssayExperiment)

# READ Gene expression
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

# READ DNA methylation
query <- GDCquery(
  project = "TCGA-READ",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value",
  sample.type = "Primary Tumor")
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

# merge gene expression and DNA methylation for READ
colnames_dataFilt_4 <- colnames(dataFilt_4)
colnames_normFeatureMatrix <- colnames(normFeatureMatrix)
participant_ids_dataFilt_4 <- substr(colnames_dataFilt_4, 9, 12)
participant_ids_normFeatureMatrix <- substr(colnames_normFeatureMatrix, 9, 12)
common_participant_ids <- intersect(participant_ids_dataFilt_4, participant_ids_normFeatureMatrix)
common_columns_dataFilt_4 <- colnames_dataFilt_4[participant_ids_dataFilt_4 %in% common_participant_ids]
common_columns_normFeatureMatrix <- colnames_normFeatureMatrix[participant_ids_normFeatureMatrix %in% common_participant_ids]
# dataFilt_4_common <- dataFilt_4[, common_columns_dataFilt_4]
# normFeatureMatrix_common <- normFeatureMatrix[, common_columns_normFeatureMatrix]
common_columns_dataFilt_4_sorted <- sort(common_columns_dataFilt_4)
common_columns_normFeatureMatrix_sorted <- sort(common_columns_normFeatureMatrix)
dataFilt_4_common_sorted <- dataFilt_4_common[, common_columns_dataFilt_4_sorted]
normFeatureMatrix_common_sorted <- normFeatureMatrix_common[, common_columns_normFeatureMatrix_sorted]
merged_data <- rbind(dataFilt_4_common_sorted, normFeatureMatrix_common_sorted)
merged_READ <- t(merged_data)
dim(merged_READ)

# COAD gene expression
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
sampleTypes_8 <- c(rep('COAD', ncol(dataFilt_8)))


# COAD DNA methylation
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value",
  sample.type = "Primary Tumor")
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

# merge gene expression and DNA methylation for COAD
colnames_dataFilt_8 <- colnames(dataFilt_8)
colnames_normFeatureMatrix_1 <- colnames(normFeatureMatrix_1)
participant_ids_dataFilt_8 <- substr(colnames_dataFilt_8, 9, 12)
participant_ids_normFeatureMatrix_1 <- substr(colnames_normFeatureMatrix_1, 9, 12)
gene_idx <- as.data.frame(t(participant_ids_dataFilt_8))
colnames(gene_idx) <- colnames_dataFilt_8
dataFilt_8_1 <- rbind(gene_idx, dataFilt_8)
met_idx <- as.data.frame(t(participant_ids_normFeatureMatrix_1))
colnames(met_idx) <- colnames_normFeatureMatrix_1
normFeatureMatrix_1_1 <-rbind(met_idx, normFeatureMatrix_1)

common_participant_ids_1 <- intersect(participant_ids_dataFilt_8, participant_ids_normFeatureMatrix_1)
common_columns_dataFilt_8 <- colnames_dataFilt_8[participant_ids_dataFilt_8 %in% common_participant_ids_1]
common_columns_normFeatureMatrix_1 <- colnames_normFeatureMatrix_1[participant_ids_normFeatureMatrix_1 %in% common_participant_ids_1]
dataFilt_8_common <- dataFilt_8_1[, common_columns_dataFilt_8]
normFeatureMatrix_common_1 <- normFeatureMatrix_1_1[, common_columns_normFeatureMatrix_1]
dataFilt_8_common <- t(dataFilt_8_common)
normFeatureMatrix_common_1 <- t(normFeatureMatrix_common_1)
merged_data_1 <- merge(dataFilt_8_common, normFeatureMatrix_common_1, by = "1")
merged_COAD <- merged_data_1[,-1]

common_cols <- Reduce(intersect, list(colnames(merged_READ), colnames(merged_COAD)))
fiREAD <- merged_READ[, common_cols]
fiCOAD <- merged_COAD[, common_cols]
X <- rbind(fiCOAD, fiREAD)
dim(X)
str(X)
X[] <- lapply(X, function(col) as.numeric(as.character(col)))
sampleTypes_COAD <- c(rep('COAD', nrow(merged_COAD)))
sampleTypes_READ <- c(rep('READ', nrow(merged_READ)))
Y <- factor(c(sampleTypes_COAD, sampleTypes_READ))
length(Y)

k = 10
sets <- lol.xval.split(X, Y, k=k, rank.low=TRUE)

library(adabag)
classifier.name <- "adaboost"
classifier.alg <- adabag::boosting
classifier.return <- 'class'

classifier.name <- "knn"
classifier.alg <- class::knn
#classifier.opts <- list(method = "knn")
classifier.return = NaN

classifier.name <- "svm"
classifier.alg <- e1071::svm
classifier.return = NaN

algs <- list(lol.project.pca, lol.project.lrlda, spp, lol.project.pls,
             lol.project.lol, spp)
names(algs) <- c("PCA", "rrLDA", 'MCDR', 'PLS', 'LOL', 'JCDR')
alg.opts=list(list(), list(), list(method='xi'), list(), list(), list(method='pca+pls'))
names(alg.opts) <- c("PCA", "rrLDA", 'MCDR', 'PLS', 'LOL', 'JCDR')


experiments <- list()
counter <- 1

X <- as.matrix(X); Y <- as.factor(Y)
n <- dim(X)[1]; d <- dim(X)[2]

len.set <- sapply(sets, function(set) length(set$train))
#max.r <- min(c(d - 1, min(len.set) - 1))
max.r <- 50
n <- length(Y)
x.n <- dim(X)[1]
if (n != x.n) {
  stop(sprintf("Your X has %d examples, but you only provide %d labels.", x.n, n))
}



log.seq <- function(from=0, to=15, length=rlen) {
  round(exp(seq(from=log(from), to=log(to), length.out=length)))
}

rs <- unique(log.seq(from=1, to=max.r, length=rlen))

results <- data.frame(exp=c(), alg=c(), xv=c(), n=c(), ntrain=c(), d=c(), K=c(), fold=c(), r=c(), lhat=c())
for (i in 1:length(algs)) {
  classifier.ret <- classifier.return
  if (classifier.name == "lda") {
    classifier.ret = "class"
    classifier.alg = MASS::lda
    if (names(algs)[i] == "QOQ") {
      classifier.alg=MASS::qda
      classifier.ret = "class"
    } else if (names(algs)[i] == "CCA") {
      classifier.alg = lol.classify.nearestCentroid
      classifier.ret = NaN
    }
  }
  tryCatch({
    xv_res <- spp.xval.optimal_dimselect(X, Y, rs, algs[[i]], sets=sets,
                                         alg.opts=alg.opts[[names(algs)[i]]], alg.return="A", classifier=classifier.alg,
                                         classifier.return=classifier.ret, k=k)
    results <- rbind(results, data.frame(sim='TCGA_READCOAD', fold=xv_res$folds.data$fold, alg=names(algs)[i],
                                         r=xv_res$folds.data$r, lhat=xv_res$folds.data$lhat))
  }, error=function(e) {print(e); return(NULL)})
}



resultso <- do.call(rbind, results)

results <- t(resultso)
write.csv(results, file=paste('/Users/ruijuanzhong/JCDR/docs/spp-paper/realdata/', 'TCGA.csv', sep=""), row.names=FALSE)
