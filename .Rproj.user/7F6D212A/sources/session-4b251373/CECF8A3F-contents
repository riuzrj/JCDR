library(datamicroarray)
library('MASS')
library('e1071')
library(parallel)
require(lolR)
library(spp)

data('sun', package = 'datamicroarray')
X = as.matrix(sun$x)
X = X[,colSums(is.na(X)) == 0]
Y = as.factor(sun$y)
K = 4


file_source <- list.files('/Users/ruijuanzhong/spp/R')
for (i in 1:length(file_source)) {
  file_source[i] <- paste('/Users/ruijuanzhong/spp/R/',file_source[i],sep ="")}
sapply(file_source,source)
file_source <- list.files('/Users/ruijuanzhong/spp/docs/lolR')
for (i in 1:length(file_source)) {
  file_source[i] <- paste('/Users/ruijuanzhong/spp/docs/lolR/',file_source[i],sep ="")}
sapply(file_source,source)

#no_cores = detectCores()
classifier.name <- "lda"
classifier.alg <- MASS::lda
classifier.return = 'class'

rlen <- 50

k = 10
sets <- lol.xval.split(X, Y, k=k, rank.low=TRUE)

experiments <- list()
counter <- 1

X <- as.matrix(X); Y <- as.factor(Y)
n <- dim(X)[1]; d <- dim(X)[2]

#pca+pls
algs <- list(spp)
names(algs) <- c('pca+pls')
algs.opts=list(list(method='pca+pls'))
names(algs.opts) <- c('pca+pls')

len.set <- sapply(sets, function(set) length(set$train))
max.r <- min(c(d - 1, min(len.set) - 1, 70))
# for pca+pls

alg.dimname="r"
alg.structured=TRUE
alg = algs[[1]]
alg.embedding="A"
alg.opts=algs.opts[[names(algs)]]
alg.return="A"
classifier=classifier.alg
classifier.opts=list()
classifier.return='class'
#k=10
# hyperparameters are the number of embedding dimensions and other options requested.
dim.embed <- list()
dim.embed[[alg.dimname]] <- max.r
alg.hparams <- c(dim.embed, alg.opts)

# tuning two parameters p and q
# p = seq(0, 24, by = 1)
# q = seq(0, 20, by = 1)

# for real data
p = seq(0, 30, by = 2)
q = seq(0, 30, by = 2)
# p = rs
# q = rs


# generate grid
grid = expand.grid(x = p, y = q)

# cross validation
cv = lapply(1:length(sets), function(i) {
  #print(i)
  set <- sets[[i]]
  #print('set')
  # if the algorithm is appropriately structured, we can avoid computing on every iteration and
  # just compute on the maximum number of embedding dimensions

  if (alg.structured) {
    # learn the projection with the algorithm specified

    mod <- do.call(alg, c(list(X=X[set$train,,drop=FALSE], Y=as.factor(Y[set$train,drop=FALSE])), alg.hparams))
    #print('dim')
    # assign the embedding  dimension
    if (is.nan(alg.embedding)) {
      A <- mod
    } else {
      A <- mod[[alg.embedding]]}
    #print(dim(A))
  }
  # calculate the error rate at each grid
  Lhat.data = apply(grid, 1, function(x){
    #print(x)
    #print(grid)
    tryCatch({
      if(x[2]==0){

        A.pq = A[,1:(x[1]+K-1)]
        #print(A.pq)
      }else{
        A.p = A[,1:(x[1]+K-1)]
        A.q = A[,(ceiling(max.r/2)+K):(ceiling(max.r/2)+K-1+x[2])]
        A.pq = cbind(A.p, A.q)}
      #print(paste('A.pq',dim(A.pq)))
      # embed the test points with the embedding matrix computed on the training data
      X.test.proj <- lol.embed(X[set$test,,drop=FALSE], A.pq)
      #print('x.test.proj')
      # project the data with the projection just learned
      # compute the trained classifier
      trained_classifier <- do.call(classifier,
                                    c(list(lol.embed(X[set$train,,drop=FALSE], A.pq),
                                           as.factor(Y[set$train,drop=FALSE])),classifier.opts))
      #print('trained.classifier')
      # and compute the held-out error for particular fold
      if (is.nan(classifier.return)) {
        Yhat <- predict(trained_classifier, X.test.proj)}
      else {Yhat <- predict(trained_classifier, X.test.proj)[[classifier.return]]}
      #print(1 - sum(Yhat == Y[set$test,drop=FALSE])/length(Yhat))
      return(data.frame(lhat=1 - sum(Yhat == Y[set$test,drop=FALSE])/length(Yhat),
                        grid_1=x[1], grid_2=x[2], r=x[1]+x[2]+K-1, fold=i))

    },
    error=function(e){print(e); return(NULL)})
  })
  #skip nulls
  Lhat.data <- Lhat.data[!sapply(Lhat.data, is.null)]
  Lhat.data = do.call(rbind, Lhat.data)
  # cv <- cv[!sapply(cv, is.null)]
  # cv = do.call(rbind, cv)
  # cv.mean = mean(cv$lhat)
  #print(cv.mean)
  return(Lhat.data)

})

results <- do.call(rbind, cv)
#write.table(Y, "output.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
results.means = aggregate(results$lhat, by = list(results$grid_1,results$grid_2), FUN = mean)
#print(dim(results.means))
min_index <- which.min(results.means$x)
#print(min_index)
min_x <- results.means$Group.1[min_index]
min_y <- results.means$Group.2[min_index]
#print('minxy')
optimal.lhat <- results.means$x[min_index]

#use optimal parameters to train the model
dim.embed <- list()
dim.embed[[alg.dimname]] <- max.r
alg.hparams <- c(dim.embed, alg.opts)
mod <- do.call(alg, c(list(X, Y), alg.hparams))
A <- mod[[alg.embedding]]
A.p = A[,1:(min_x+K-1)]
A.pq = cbind(A[, 1:(min_x+K-1)], A[, (ceiling(max.r/2)+K):(ceiling(max.r/2)+K-1+min_y)])

trained_classifier <- do.call(classifier,
                              c(list(lol.embed(X, A.pq),
                                     as.factor(Y)),classifier.opts))



bootstrap_sample <- function(data) {
  sampled_rows <- sample(nrow(data), replace = TRUE)
  data[sampled_rows, ]
}

generate_bootstrap_samples <- function(data, n_samples) {
  lapply(1:n_samples, function(i) bootstrap_sample(data))
}

bootstrap_samples <- generate_bootstrap_samples(X, 100)

iteration_function <- function(j) {

  X_sim <- bootstrap_samples[[j]]
  #add noise
  # noise_range <- 0.1
  # X_sim <- X_sim + runif(length(X_sim), min = -noise_range, max = noise_range)
  noise.level = 0.1
  noise = matrix(rnorm(nrow(X_sim) * ncol(X_sim), sd = noise.level), nrow(X_sim), ncol(X_sim))
  X_sim = X_sim + noise
  generated_Y <- predict(trained_classifier, lol.embed(X_sim, A.pq))[[classifier.return]]

  Y = generated_Y
  X = X_sim

  k = 10
  sets <- lol.xval.split(X, Y, k=k, rank.low=TRUE)

  algs <- list(lol.project.pca, lol.project.lrlda, lol.project.pls,
               lol.project.lol, spp, spp)
  algs.name <- c("PCA", "rrLDA", "PLS", "LOL", 'pca+pls', 'Naive_pca+pls')
  names(algs) <-algs.name
  algs.opts=list(list(), list(), list(),
                 list(), list(method='pca+pls'), list(method='xi'))
  names(algs.opts) <- algs.name

  # algs <- list(lol.project.lol)
  # algs.name <- c('pca+pls')
  # names(algs) <- algs.name
  # algs.opts <- list(list(method='pca+pls'))
  # names(algs.opts) <- algs.name

  experiments <- list()
  counter <- 1

  X <- as.matrix(X); Y <- as.factor(Y)
  n <- dim(X)[1]; d <- dim(X)[2]

  len.set <- sapply(sets, function(set) length(set$train))
  max.r <- min(c(d - 1, min(len.set) - 2, 60))

  n <- length(Y)
  x.n <- dim(X)[1]
  if (n != x.n) {
    stop(sprintf("Your X has %d examples, but you only provide %d labels.", x.n, n))
  }

  log.seq <- function(from=0, to=15, length=rlen) {
    round(exp(seq(from=log(from), to=log(to), length.out=length)))
  }

  rs <- unique(log.seq(from=1, to=max.r, length=rlen))

  results <- data.frame(sim=c(), iter=c(),fold=c(),alg=c(),r=c(),lhat=c())
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
      }else if (names(algs)[i]=='xi'){
        classifier.alg = e1071::svm
        classifier.ret = "class"
      }
    }
    tryCatch({
      xv_res <- spp.xval.optimal_dimselect(X, Y, rs, algs[[i]], sets=sets,
                                           alg.opts=algs.opts[[names(algs)[i]]], alg.return="A", classifier=classifier.alg,
                                           classifier.return=classifier.ret, k=k)
      results <- rbind(results, data.frame(sim='uni_nois', iter=j,
                                           fold=xv_res$folds.data$fold,
                                           alg=names(algs)[i],
                                           r=xv_res$folds.data$r,
                                           lhat=xv_res$folds.data$lhat))
      results = unique(results)
      #print(results)
    }, error=function(e) {print(e); return(NULL)})
  }

  return(results)
  #write.csv(results, file = sprintf('sim_%d', j))
}


results_all <- mclapply(1:100, iteration_function, mc.cores = 100)
#result_all <- lapply(1:1, iteration_function)

results_combined <- do.call(rbind, results_all)

write.csv(results_combined,'gaussian_noise_sun.csv')

results_means = mean(results_combined[,2])
