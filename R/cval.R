#' Embedding Cross Validation
#'
#' A function for performing cross-validation for a given emdedding model.
#' This function produces fold-wise cross-validated misclassification rates.
#' Users can optionally specify custom embedding techniques with proper configuration of \code{alg.*} parameters and hyperparameters.
#'
#' @importFrom MASS lda
#' @importFrom stats predict
#' @importFrom lolR lol.xval.split
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the number of embedding dimensions desired, where \code{r <= d}.
#' @param alg the algorithm to use for embedding. Should be a function that accepts inputs \code{X}, \code{Y}, and has a parameter for \code{alg.dimname} if \code{alg} is supervised, or just \code{X} and \code{alg.dimname} if \code{alg} is unsupervised.This algorithm should return a list containing a matrix that embeds from {d} to {r <= d} dimensions.
#' @param sets a user-defined cross-validation set. Defaults to \code{NULL}.
#' \itemize{
#' \item{\code{is.null(sets)} randomly partition the inputs \code{X} and \code{Y} into training and testing sets.}
#' \item{\code{!is.null(sets)} use a user-defined partitioning of the inputs \code{X} and \code{Y} into training and testing sets. Should be in the format of the outputs from \code{\link{lol.xval.split}}. That is, a \code{list} with each element containing \code{X.train}, an \code{[n-k][d]} subset of data to test on, \code{Y.train}, an \code{[n-k]} subset of class labels for \code{X.train}; \code{X.test}, an \code{[n-k][d]} subset of data to test the model on, \code{Y.train}, an \code{[k]} subset of class labels for \code{X.test}.}
#' }
#' @param alg.dimname the name of the parameter accepted by \code{alg} for indicating the embedding dimensionality desired. Defaults to \code{r}.
#' @param alg.opts the hyper-parameter options you want to pass into your algorithm, as a keyworded list. Defaults to \code{list()}, or no hyper-parameters.
#' @param alg.embedding the attribute returned by \code{alg} containing the embedding matrix. Defaults to assuming that \code{alg} returns an embgedding matrix as \code{"A"}.
#' \itemize{
#' \item \code{!is.nan(alg.embedding)} Assumes that \code{alg} will return a list containing an attribute, \code{alg.embedding}, a \code{[d, r]} matrix that embeds \code{[n, d]} data from \code{[d]} to \code{[r < d]} dimensions.
#' \item \code{is.nan(alg.embedding)} Assumes that \code{alg} returns a \code{[d, r]} matrix that embeds \code{[n, d]} data from \code{[d]} to \code{[r < d]} dimensions.
#' }
#'@param classifier the classifier to use for assessing performance. The classifier should accept \code{X}, a \code{[n, d]} array as the first input, and \code{Y}, a \code{[n]} array of labels, as the first 2 arguments. The class should implement a predict function, \code{predict.classifier}, that is compatible with the \code{stats::predict} \code{S3} method. Defaults to \code{MASS::lda}.
#' @param classifier.opts any extraneous options to be passed to the classifier function, as a list. Defaults to an empty list.
#' @param classifier.return if the return type is a list, \code{class} encodes the attribute containing the prediction labels from \code{stats::predict}. Defaults to the return type of \code{MASS::lda}, \code{class}.
#' \itemize{
#' \item{\code{!is.nan(classifier.return)} Assumes that \code{predict.classifier} will return a list containing an attribute, \code{classifier.return}, that encodes the predicted labels.}
#' \item{\code{is.nan(classifier.return)} Assumes that \code{predict.classifer} returns a \code{[n]} vector/array containing the prediction labels for \code{[n, d]} inputs.}
#' }
#' @param k the cross-validated method to perform. Defaults to \code{'loo'}. If \code{sets} is provided, this option is ignored. See \code{\link{lol.xval.split}} for details.
#' \itemize{
#' \item{\code{'loo'} Leave-one-out cross validation}
#' \item{\code{isinteger(k)}  perform \code{k}-fold cross-validation with \code{k} as the number of folds.}
#' }
#' @param rank.low whether to force the training set to low-rank. Defaults to \code{FALSE}. If \code{sets} is provided, this option is ignored. See \code{\link{lol.xval.split}} for details.
#' \itemize{
#' \item{if \code{rank.low == FALSE}, uses default cross-validation method with standard \code{k}-fold validation. Training sets are \code{k-1} folds, and testing sets are \code{1} fold, where the fold held-out for testing is rotated to ensure no dependence of potential downstream inference in the cross-validated misclassification rates.}
#' \item{if \code{rank.low == TRUE}, users cross-validation method with \code{ntrain = min((k-1)/k*n, d)} sample training sets, where \code{d}  is the number of dimensions in \code{X}. This ensures that the training data is always low-rank, \code{ntrain < d + 1}. Note that the resulting training sets may have \code{ntrain < (k-1)/k*n}, but the resulting testing sets will always be properly rotated \code{ntest = n/k} to ensure no dependencies in fold-wise testing.}
#' }
#' @param ... additional arguments to be passed to the \code{alg} function.
#' @return a list containing the following elements:
#' \itemize{
#' \item{\code{folds.data}}{the results of the per-fold accuracy of the classifier.}
#' \item{\code{foldmeans.data}}{the mean accuracy of the classifier across all folds.}
#' \item{\code{optimal.lhat}}{the optimal classifier accuracy.}
#' \item{\code{optimal.r}}{the optimal embedding dimensionality.}
#' \item{\code{model}}{the model object returned by the classifier at the optimal embedding dimensionality.}
#'
#' @examples
#' library(spp)
#' library(MASS)
#' library(lolR)
#' data(iris)
#' X <- iris[,1:4]; Y <- iris[,5]
#' cval <- spp.xval.optimal_dimselect(X, Y, rs=c(1,2), spp, alg.opts=list(method='mvi'),
#'                                    alg.embedding="A", classifier.opts=list(),k=10)
#'
spp.xval.optimal_dimselect <- function(X, Y, rs, alg, sets=NULL, alg.dimname="r", alg.opts=list(), alg.embedding="A",
                                       alg.structured=TRUE, classifier=lda, classifier.opts=list(),
                                       classifier.return="class", k='loo', rank.low=FALSE, ...) {
  if (length(alg.opts) == 0){

    #write.table(X, "output.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
    d <- dim(X)[2]
    Y <- factor(Y)
    n <- length(Y)
    x.n <- dim(X)[1]
    if (n != x.n) {
      stop(sprintf("Your X has %d examples, but you only provide %d labels.", x.n, n))
    }
    # check that if the user specifies the cross-validation set, if so, that it is correctly set up
    # otherwise, do it for them
    if (is.null(sets)) {
      sets <- lol.xval.split(X, Y, k=k, rank.low=rank.low)
    } else {
      lol.xval.check_xv_set(sets, n)
    }

    # compute the top r embedding dimensions
    max.r <- max(rs)

    # hyperparameters are the number of embedding dimensions and other options requested.
    dim.embed <- list()
    dim.embed[[alg.dimname]] <- max.r
    alg.hparams <- c(dim.embed, alg.opts)


    Lhat.data <- lapply(1:length(sets), function(i) {
      set <- sets[[i]]
      # if the algorithm is appropriately structured, we can avoid computing on every iteration and
      # just compute on the maximum number of embedding dimensions

      if (alg.structured) {
        # learn the projection with the algorithm specified
        mod <- do.call(alg, c(list(X=X[set$train,,drop=FALSE], Y=as.factor(Y[set$train,drop=FALSE])), alg.hparams))
        # assign the embedding  dimension
        if (is.nan(alg.embedding)) {
          A <- mod
        } else {
          A <- mod[[alg.embedding]]}
      }


      # check fold with every desired possibility of embedding dimensions
      res.rs <- lapply(rs, function(r){
        tryCatch({
          # if appropriately structured just take the top r dimensions of the embedding computed on max.r
          if (alg.structured) {
            A.r <- A[, 1:r]}
          else {
            # otherwise, compute A.r on the new embedding dimension every time
            alg.hparams[[alg.dimname]] <- r
            mod <- do.call(alg, c(list(X=X[set$train,,drop=FALSE], Y=as.factor(Y[set$train,drop=FALSE])),alg.hparams))
            if (is.nan(alg.embedding)) {
              A.r <- mod}
            else {
              A.r <- mod[[alg.embedding]]}}
          # embed the test points with the embedding matrix computed on the training data
          X.test.proj <- lol.embed(X[set$test,,drop=FALSE], A.r)
          # project the data with the projection just learned
          # compute the trained classifier
          trained_classifier <- do.call(classifier,
                                        c(list(lol.embed(X[set$train,,drop=FALSE], A.r),
                                               as.factor(Y[set$train,drop=FALSE])),classifier.opts))

          # and compute the held-out error for particular fold
          if (is.nan(classifier.return)) {
            Yhat <- predict(trained_classifier, X.test.proj)}
          else {Yhat <- predict(trained_classifier, X.test.proj)[[classifier.return]]}

          return(data.frame(lhat=1 - sum(Yhat == Y[set$test,drop=FALSE])/length(Yhat), r=r, fold=i, Yhat))},
          error=function(e){print(e); return(NULL)})
      })
      # skip nulls
      res.rs <- res.rs[!sapply(res.rs, is.null)]
      res.rs = do.call(rbind, res.rs)
      return(res.rs)
    })
    results <- do.call(rbind, Lhat.data)

    # compute fold-wise average
    results.means <- aggregate(lhat ~ r, data = results, FUN = nan.mean)

    # find the number of embedding dimensions with lowest average error
    optimal.idx <- which(results.means$lhat == min(results.means$lhat))
    best.r <- min(results.means$r[optimal.idx])  # best is the minimum of the possible choices
    # recompute on all data with the best number of embedding dimensions
    alg.hparams[[alg.dimname]] <- best.r
    model <- do.call(alg, c(list(X=X, Y=Y), alg.hparams))

    if (is.nan(alg.embedding)) {
      A <- model}
    else {
      A <- model[[alg.embedding]]}

    # and train the classifier for good measure
    class <- do.call(classifier, c(list(lol.embed(X, A), Y), classifier.opts))
    return(list(folds.data=results,
                foldmeans.data=results.means,
                optimal.lhat=results.means$lhat[optimal.idx],
                optimal.r=best.r,
                model=model, classifier=class))
  }
   if (alg.opts == 'mvi'|| alg.opts == 'xi'){

    # print("MVI")
    d <- dim(X)[2]
    Y <- factor(Y)
    n <- length(Y)
    x.n <- dim(X)[1]
    if (n != x.n) {
      stop(sprintf("Your X has %d examples, but you only provide %d labels.", x.n, n))
    }
    # check that if the user specifies the cross-validation set, if so, that it is correctly set up
    # otherwise, do it for them
    if (is.null(sets)) {
      sets <- lolR::lol.xval.split(X, Y, k=k, rank.low=rank.low)
    } else {
      lolR::lol.xval.check_xv_set(sets, n)
    }

    # compute the top r embedding dimensions
    max.r <- max(rs)

    # hyperparameters are the number of embedding dimensions and other options requested.
    dim.embed <- list()
    dim.embed[[alg.dimname]] <- max.r
    alg.hparams <- c(dim.embed, alg.opts)


    Lhat.data <- lapply(1:length(sets), function(i) {
      #print(i)
      set <- sets[[i]]
      # if the algorithm is appropriately structured, we can avoid computing on every iteration and
      # just compute on the maximum number of embedding dimensions

      if (alg.structured) {
        # learn the projection with the algorithm specified
        mod <- do.call(alg, c(list(X=X[set$train,,drop=FALSE], Y=as.factor(Y[set$train,drop=FALSE])), alg.hparams))
        # assign the embedding  dimension
        if (is.nan(alg.embedding)) {
          A <- mod
        } else {
          A <- mod[[alg.embedding]]}
      }

      # check fold with every desired possibility of embedding dimensions
      res.rs <- lapply(rs, function(r){
        tryCatch({
          # if appropriately structured just take the top r dimensions of the embedding computed on max.r
          if (alg.structured) {
            A.r <- A[, 1:r]}
          else {
            # otherwise, compute A.r on the new embedding dimension every time
            alg.hparams[[alg.dimname]] <- r
            mod <- do.call(alg, c(list(X=X[set$train,,drop=FALSE], Y=as.factor(Y[set$train,drop=FALSE])),alg.hparams))
            if (is.nan(alg.embedding)) {
              A.r <- mod}
            else {
              A.r <- mod[[alg.embedding]]}}
          # embed the test points with the embedding matrix computed on the training data
          X.test.proj <- lol.embed(X[set$test,,drop=FALSE], A.r)
          # project the data with the projection just learned
          # compute the trained classifier
          trained_classifier <- do.call(classifier,
                                        c(list(lol.embed(X[set$train,,drop=FALSE], A.r),
                                               as.factor(Y[set$train,drop=FALSE])),classifier.opts))
          #print('5')
          # and compute the held-out error for particular fold
          if (is.nan(classifier.return)) {
            Yhat <- predict(trained_classifier, X.test.proj)}
          else {Yhat <- predict(trained_classifier, X.test.proj)[[classifier.return]]}
          return(data.frame(lhat=1 - sum(Yhat == Y[set$test,drop=FALSE])/length(Yhat), r=r, fold=i))},
          error=function(e){print(e); return(NULL)})
      })
      # skip nulls
      res.rs <- res.rs[!sapply(res.rs, is.null)]
      res.rs = do.call(rbind, res.rs)
      #print('6')
      return(res.rs)
    })
    results <- do.call(rbind, Lhat.data)

    # compute fold-wise average
    results.means <- aggregate(lhat ~ r, data = results, FUN = nan.mean)

    # find the number of embedding dimensions with lowest average error
    optimal.idx <- which(results.means$lhat == min(results.means$lhat))
    best.r <- min(results.means$r[optimal.idx])  # best is the minimum of the possible choices
    # recompute on all data with the best number of embedding dimensions
    alg.hparams[[alg.dimname]] <- best.r
    model <- do.call(alg, c(list(X=X, Y=Y), alg.hparams))

    if (is.nan(alg.embedding)) {
      A <- model}
    else {
      A <- model[[alg.embedding]]}

    # and train the classifier for good measure
    class <- do.call(classifier, c(list(lol.embed(X, A), Y), classifier.opts))
    #print('MVI2')
    return(list(folds.data=results,
                foldmeans.data=results.means,
                optimal.lhat=results.means$lhat[optimal.idx],
                optimal.r=best.r,
                model=model, classifier=class))
  }


  ## pca+pls part
  if (alg.opts=='pca+pls'){

    #print('pca+pls')
    info <- lol.utils.info(X, Y, robust=FALSE)
    #print('pca+pls')
    K <- info$K
    d <- dim(X)[2]
    Y <- factor(Y)
    n <- length(Y)
    x.n <- dim(X)[1]
    if (n != x.n) {
      stop(sprintf("Your X has %d examples, but you only provide %d labels.", x.n, n))
    }
    # check that if the user specifies the cross-validation set, if so, that it is correctly set up
    # otherwise, do it for them
    if (is.null(sets)) {
      sets <- lol.xval.split(X, Y, k=k, rank.low=rank.low)
    } else {
      lol.xval.check_xv_set(sets, n)
    }

    # compute the top r embedding dimensions
    max.r <- max(rs)+30

    # hyperparameters are the number of embedding dimensions and other options requested.
    dim.embed <- list()
    dim.embed[[alg.dimname]] <- max.r
    alg.hparams <- c(dim.embed, alg.opts)

    # tuning two parameters p and q
    # p = seq(0, 24, by = 1)
    # q = seq(0, 20, by = 1)

    # for real data
    # p = seq(0, 35, by = 1)
    # q = seq(0, 35, by = 1)

    p = seq(1, 30, by = 1)
    q = seq(0, 20, by = 1)
    # generate grid
    grid = expand.grid(x = p, y = q)

    # calculate the error rate at each grid


    # cross valisation
    cv = lapply(1:length(sets), function(i) {
      #print(i)
      set <- sets[[i]]
      # if the algorithm is appropriately structured, we can avoid computing on every iteration and
      # just compute on the maximum number of embedding dimensions

      if (alg.structured) {
        # learn the projection with the algorithm specified
        mod <- do.call(alg, c(list(X=X[set$train,,drop=FALSE], Y=as.factor(Y[set$train,drop=FALSE])), alg.hparams))
        # assign the embedding  dimension
        if (is.nan(alg.embedding)) {
          A <- mod
        } else {
          A <- mod[[alg.embedding]]}
      }

      Lhat.data = apply(grid, 1, function(x){
        #print('grid')
        tryCatch({
          if(x[2]==0){
            A.pq = A[,1:(x[1]+K-1)]
            #A.p = A[,1:x[1]]
          }else{
            A.p = A[,1:(x[1]+K-1)]
            #A.p = A[,1:x[1]]
            A.q = A[,(ceiling(max.r/2)+K):(ceiling(max.r/2)+K-1+x[2])]
            A.pq = cbind(A.p, A.q)

            X_proj <- lol.embed(X[set$train,,drop=FALSE], A.pq)

            screen_result <- VariableScreening::screenIID(X_proj,as.factor(Y[set$train, drop=FALSE]),method = 'MV-SIS')
            idx = screen_result$rank
            A.pq = A.pq[,idx]
          }
          # embed the test points with the embedding matrix computed on the training data
          X.test.proj <- lol.embed(X[set$test,,drop=FALSE], A.pq)
          # project the data with the projection just learned
          # compute the trained classifier
          trained_classifier <- do.call(classifier,
                                        c(list(lol.embed(X[set$train,,drop=FALSE], A.pq),
                                               as.factor(Y[set$train,drop=FALSE])),classifier.opts))

          # and compute the held-out error for particular fold
          if (is.nan(classifier.return)) {
            Yhat <- predict(trained_classifier, X.test.proj)}
          else {Yhat <- predict(trained_classifier, X.test.proj)[[classifier.return]]}

          return(data.frame(lhat=1 - sum(Yhat == Y[set$test,drop=FALSE])/length(Yhat),
                            grid_1=x[1],grid_2=x[2],
                            #r=x[1]+x[2],
                            r=x[1]+x[2]+K-1,
                            fold=i, Yhat=Yhat, Y=Y[set$test,drop=FALSE]))
          # return(data.frame(lhat=1 - sum(Yhat == Y[set$test,drop=FALSE])/length(Yhat),
          #                   grid_1=x[1],grid_2=x[2],
          #                   #r=x[1]+x[2],
          #                   r=x[1]+x[2]+K-1,
          #                   fold=i))
        },
        error=function(e){print(e); return(NULL)})
      })
      #skip nulls
      Lhat.data <- Lhat.data[!sapply(Lhat.data, is.null)]
      Lhat.data = do.call(rbind, Lhat.data)
      return(Lhat.data)
    })

    results <- do.call(rbind, cv)
    #print(results[1,])
    results.means = aggregate(results$lhat, by = list(results$grid_1,results$grid_2), FUN = mean)
    #print(results.means[1,])
    results.means$r = results.means$Group.1 + results.means$Group.2 + K - 1
    names(results.means)[1:3] <- c("Group.1", "Group.2", "lhat")

    results_df <- results.means %>%
      group_by(r) %>%
      filter(lhat == min(lhat)) %>%
      ungroup() %>%
      dplyr::select(Group.1, Group.2) %>%
      distinct()

    results <- results %>%
      semi_join(results_df, by = c("grid_1" = "Group.1", "grid_2" = "Group.2"))

    #print(results[1,])

    #print(dim(results.means))
    min_index <- which.min(results.means$x)
    #print(min_index)
    min_x <- results.means$Group.1[min_index]
    min_y <- results.means$Group.2[min_index]
    #print('minxy')
    optimal.lhat <- results.means$x[min_index]

    ##  best.r??
    #print(optimal.lhat)
    ### need to correct
    model <- do.call(alg, c(list(X=X, Y=Y), alg.hparams))

    if (is.nan(alg.embedding)) {
      A <- model}
    else {
      A <- model[[alg.embedding]]}

    # and train the classifier for good measure
    class <- do.call(classifier, c(list(lol.embed(X, A), Y), classifier.opts))

    return(list(folds.data=results,
                foldmeans=results.means,
                optimal.lhat=optimal.lhat,
                optimal.r=paste(min_x+min_y+1),
                model=model, classifier=class))
  }
}

