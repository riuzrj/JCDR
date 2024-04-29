# Parallelize Stuff
#=========================#
library(dimRed)
require(MASS)
library(parallel)
library(datamicroarray)
library(dplyr)

file_source <- list.files('/Users/ruijuanzhong/spp/docs/lolR')
for (i in 1:length(file_source)) {
  file_source[i] <- paste('/Users/ruijuanzhong/spp/docs/lolR/',file_source[i],sep ="")}
sapply(file_source,source)

no_cores = detectCores()
classifier.name <- "lda"
classifier.alg <- MASS::lda
classifier.return = 'class'

rlen <- 50

# Setup Algorithms
#==========================#
# algs <- list(lol.project.pca, lol.project.lrlda, lol.project.lrcca, lol.project.rp, lol.project.pls,
#              lol.project.lol)
# names(algs) <- c("PCA", "LRLDA", "CCA", "RP", "PLS", "LOL")
# alg.opts=list(list(), list(), list(), list(), list(), list(), list(second.moment="quadratic"))
# names(alg.opts) <- c("PCA", "LRLDA", "CCA", "RP", "PLS", "LOL", "QOL")

#data load
#data('yeoh', package = 'datamicroarray')
#data('burczynski', package = 'datamicroarray')
#data('gravier', package = 'datamicroarray')
data('nakayama', package = 'datamicroarray')
#data('chowdary', package = 'datamicroarray')
#data('chin', package = 'datamicroarray')
#data('gordon', package = 'datamicroarray')
#data('singh', package = 'datamicroarray')
#data('alon', package = 'datamicroarray')
#data('west', package = 'datamicroarray')
#data('tian', package = 'datamicroarray')
#data('sun', package = 'datamicroarray')

# data = read.csv('/Users/ruijuanzhong/Project/Real_data/Breast_GSE45827.csv')
# data = data[,-1]
# X = data[,2:ncol(data)]
# Y = data[,1]
# data = t(data)
# data = data[-1,]
# X = data[,2:ncol(data)]
# X <- matrix(as.numeric(X), nrow = nrow(X), ncol = ncol(X))
# dim(X)
# Y = data[,1]

# Y = sun$y
# X = sun$x
# #去掉缺失值所在列
# X = X[,colSums(is.na(X)) == 0]

# Y = tian$y
# X = tian$x

# Y = west$y
# X = west$x

# Y = alon$y
# X = alon$x

# Y = singh$y
# X = singh$x

# Y = gordon$y
# X = gordon$x

# Y = chin$y
# X = chin$x

# Y = chowdary$y
# X = chowdary$x

# Y = gravier$y
# X = gravier$x

# Y = burczynski$y
# X = burczynski$x
table(Y)
Y = nakayama$y
X = nakayama$x
class_counts <- table(Y)

# 找出样本量为 3 的类
classes_to_remove <- names(class_counts[class_counts == 3])

# 删除这些类
X <- X[!(Y %in% classes_to_remove), ]
Y <- Y[!(Y %in% classes_to_remove)]
Y <- factor(Y)

# Y = yeoh$y
# X = yeoh$x

# category_counts <- table(Y)
# total_count <- sum(category_counts)
# category_proportions <- category_counts / total_count
# print(category_proportions)

#remove samples with proportion < 0.1
# 从数据集中删除 "BCR" 和 "MLL" 类别的样本
# x_filtered <- yeoh$x[!(yeoh$y %in% c("BCR", "MLL")), ]
# y_filtered <- yeoh$y[!(yeoh$y %in% c("BCR", "MLL"))]
# X = as.matrix(x_filtered)
# Y = as.factor(y_filtered)

# x_filtered <- sun$x[!(sun$y %in% c("astrocytomas")), ]
# y_filtered <- sun$y[!(sun$y %in% c("astrocytomas"))]
# X = as.matrix(x_filtered)
# X = X[,colSums(is.na(X)) == 0]
# Y = as.factor(y_filtered)


#k = 'loo'  # number of folds
k = 10
sets <- lol.xval.split(X, Y, k=k, rank.low=TRUE)


algs <- list(lol.project.pca, lol.project.lrlda, spp, lol.project.pls,
             lol.project.lol, spp)
names(algs) <- c("PCA", "rrLDA", 'Naive_pca-pls', 'PLS', 'LOL', 'pca+pls')
alg.opts=list(list(), list(), list(method='xi'), list(), list(), list(method='pca+pls'))
names(alg.opts) <- c("PCA", "rrLDA", 'Naive_pca-pls', 'PLS', 'LOL', 'pca+pls')

# algs <- list(lol.project.lol)
# names(algs) <- c("pca+pls")
# alg.opts=list(list(method='pca+pls'))
# names(alg.opts) <- c("pca+pls")

# algs <- list(lol.project.lol)
# names(algs) <- c("xi")
# alg.opts=list(list(method= 'mvi'))
# names(alg.opts) <- c("xi")
#
# algs <- list(lol.project.pls, lol.project.lol)
# names(algs) <- c("PLS", "LOL")
# alg.opts=list(list(), list())
# names(alg.opts) <- c("PLS", "LOL")

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

# opath <- '/Users/ruijuanzhong/Project/lol.pro/data/singh_data/'
# dir.create(opath)

log.seq <- function(from=0, to=15, length=rlen) {
  round(exp(seq(from=log(from), to=log(to), length.out=length)))
}

rs <- unique(log.seq(from=1, to=max.r, length=rlen))
#print(length(rs))
#print(rs)
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
    results <- rbind(results, data.frame(sim='nakayama', fold=xv_res$folds.data$fold, alg=names(algs)[i],
                                         r=xv_res$folds.data$r, lhat=xv_res$folds.data$lhat))
  }, error=function(e) {print(e); return(NULL)})
}

# classifier <- "RandomGuess"
# model <- do.call(classifier.algs[[classifier]], list(X[sets[[1]]$train, ], factor(Y[sets[[1]]$train], levels=unique(Y[sets[[1]]$train]))))
# results <- rbind(results, data.frame(exp=taskname, alg=algs, xv=k, n=n, ntrain=length(sets[[1]]$train), d=d, K=length(unique(Y)),
#                                      fold=fold$fold, r=NaN, lhat=1 - max(model$priors), repo=dat$repo))
#
# results <- results[complete.cases(results$lhat),]
#
# saveRDS(results, file=paste(opath, taskname, '_', fold$fold, '.rds', sep=""))
#return(results)

#time.after=Sys.time()

resultso <- do.call(rbind, results)

results <- t(resultso)
write.csv(results, file=paste('/Users/ruijuanzhong/Project/', 'sun_aug.csv', sep=""), row.names=FALSE)
# filter out bad rows
#resultso <- resultso[complete.cases(resultso$lhat) & !(is.infinite(resultso$lhat)) & complete.cases(resultso),]
saveRDS(resultso, file.path(opath, paste(classifier.name, 'real_nakayama_nc.rds', sep="")))


library(tidyverse)
head(results)
library(dplyr)
library(tidyr) # 为了使用 unnest 函数

# 假设 desired_r 是已经定义好的
desired_r <- seq(min(results$r), 50, by = 1)

# 对数据进行分组，并应用 spline 函数进行平滑处理
filterd_results <- results %>%
  group_by(sim, alg, r) %>%
  summarize(lhat = mean(lhat)) %>%
  #filter(!(alg %in% c('pca+pls') & !(r %in% rs))) %>%
  #把大于 30 的 r 值过滤掉
  filter(r <= 50)
#   mutate(smoothed = list(spline(rs, lhat, xout = desired_r, method = 'fmm', ties = mean))) %>%
# # 使用 unnest 获得平滑后的数据
#   unnest(c(smoothed))

#用 spline 函数进行平滑处理
# smoothed_results <- filterd_results %>%
#   group_by(sim, alg) %>%
#   summarize(x = list(desired_r), y = list(spline(r, lhat, xout = desired_r,
#                                                  method = 'fmm', ties = mean)$y)) %>%
#   unnest(c(x, y))

# dt_plot <- smoothed_results %>%
#   group_by(sim, alg) %>%
#   mutate(best.lhat = min(y),
#          best.r = list(x[y == best.lhat]),
#          lhat.thresh = best.lhat * 1.05,
#          lhat.beats = y <= lhat.thresh,
#          r.star = best.r)

#用 loess 函数进行平滑处理
smoothed_results <- filterd_results %>%
  group_by(sim, alg) %>%
  do({
    # 对每个组应用 loess 平滑
    smooth_span = 0.75  # 这个值可以调整来控制平滑程度，数值范围是 0 到 1
    data.frame(r = .$r, lhat = predict(loess(lhat ~ r, data = ., span = smooth_span)))
  }) %>%
  ungroup()

dt_plot <- smoothed_results %>%
  group_by(sim, alg) %>%
  mutate(best.lhat = min(lhat),
         best.r = r[lhat == best.lhat],
         lhat.thresh = best.lhat * 1.05,
         lhat.beats = lhat <= lhat.thresh,
         #r.star = min(r[lhat.beats]))
         r.star = best.r) #%>%
#把 LRLDA 改为 rrLDA
# mutate(sim = recode_factor(sim, "uni_nois" = "gaussian_noise"))


#用 spline 函数进行平滑处理
# smoothed_results <- filterd_results %>%
#   group_by(sim, alg) %>%
#   summarize(x = list(desired_r), y = list(spline(r, lhat, xout = desired_r, method = 'fmm', ties = mean)$y)) %>%
#   unnest(c(x, y))
#
# dt_plot <- smoothed_results %>%
#   group_by(sim, alg) %>%
#   mutate(best.lhat = min(y),
#          best.r = list(x[y == best.lhat]),
#          lhat.thresh = best.lhat * 1.05,
#          lhat.beats = y <= lhat.thresh,
#          r.star = min(x[lhat.beats]))



# Continue with your plotting code
# ...

# dt_plot <- results %>%
#   group_by(sim,alg,r) %>%
#   summarise(lhat=mean(lhat)) %>%
#   #filter(!(alg %in% c('pca+pls') & !(r %in% rs))) %>%
#   mutate(best.lhat=min(lhat),
#          best.r=list(r[lhat == best.lhat]),
#          lhat.thresh=best.lhat*1.05,
#          lhat.beats=lhat <= best.lhat,
#          r.star=min(r[lhat.beats]))

#mutate(lhat_adjust = spline(rs, lhat, xout = desired_r, method = 'fmm', ties = mean))

#把 alg 中‘pca+pls’这个类中 r 大于 50 的数据去掉
#filter(!(alg %in% c('pca+pls') & r > 50))
#筛选出alg 中‘pca+pls’这个类中 r等于 rs 的数值



# mutate(alg=factor(recode_factor(alg, "QOL"="QOQ", "LRLDA"="rrLDA"),
#                   levels=algs.name,
#                   ordered=TRUE)) %>%
#filter(!((alg %in% c("RLOL")) & (sim == "Cross"))) %>%
# filter(!(alg %in% c("MVI"))) %>%
#filter(!((alg %in% c("QOQ")) & (sim %in% c("Trunk-3", "Robust")))) %>%
#mutate(
# sim=recode_factor(sim, "Trunk-3"="(A) Trunk-3", "Robust"="(B) Robust", "Cross"="(C) Cross",
# "class.diff-3"="(D) class.diff-3"),
# sim=factor(sim, levels=c("(A) Trunk-3", "(B) Robust", "(C) Cross","(D) class.diff"), ordered=TRUE))
# sim=factor(sim, levels=c("Trunk-3", "Robust", "Cross","class.diff-3"), ordered=TRUE))
#sim=factor(sim, levels=sims.name, ordered=TRUE))
head(dt_plot)
# nan.idx<- which(rowSums(is.na(dt_plot))==TRUE)


unique(dt_plot$alg)
color_lib <- c("PCA"='#1f77b4',"LOL"='#ff7f0e', "PLS"='#2ca02c',
               "rrLDA"= '#9467bd', 'pca+pls'='#d62728', 'Naive_pca-pls'= '#8c564b')


line_lib <- c("PCA"="solid", "LOL"="solid", "PLS"="solid",
              "rrLDA"="solid", 'pca+pls'='solid', 'Naive_pca-pls'='solid')

dt.star <- dt_plot %>% filter(r == r.star)
dt_plot %>%
  # slice(nan.idx) %>%
  ggplot(aes(x=r, y=lhat, color=alg)) +
  geom_line(aes(linetype=alg), linewidth=0.7) +
  scale_linetype_manual(values=line_lib, name="Algorithm") +
  scale_color_manual(values=color_lib, name="Algorithm") +
  geom_point(data=dt.star, size=3, shape=16) +
  # facet_grid(.~sim) +
  facet_wrap(.~sim,scales = 'free',nrow = 2) +
  xlab("r") +
  ylab("error") +
  guides(color=guide_legend(ncol=1), linetype=guide_legend(ncol=1))
# labs(
#   title = paste('Classifier:',classifier.name,', n=',n,', d=',dim(X[2]),', rep=')  ) +
# theme(text=element_text(size=18))
# strip.text.x=element_blank(),
# strip.background=element_blank())


