#' Supervised PCA-PLS
#'
#' A function that performs supervised PCA-PLS.
#'
#' @importFrom lolR pls
#' @importFrom lolR lrlda
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection. Note that \code{r >= K}, and \code{r < d}.
#' @param first.moment the function to capture the first moment. Defaults to \code{'delta'}.
#' \itemize{
#' \item{\code{'delta'} capture the first moment with the hyperplane separating the per-class means.}
#' \item{\code{FALSE} do not capture the first moment.}
#' }
#' @param second.moment the function to capture the second moment. Defaults to \code{'linear'}.
#' @param orthogonalize whether to orthogonalize the projection matrix. Defaults to \code{FALSE}.
#' @param second.moment.xfm whether to use extraneous options in estimation of the second moment component.
#' @param second.moment.xfm.opts optional arguments to pass to the \code{second.moment.xfm} option specified.
#' @param method the method to use for the projection. Defaults to \code{'pca+pls'}. Options include: \code{'mvi'}, \code{'xi'}, \code{'double.proj'}, and \code{'pca+pls'}.
#' @param orth.first whether to orthogonalize the projection matrix with first moment and second
#' moment components. Defaults to \code{FALSE}.
#' @param robust.first whether to use robust PCA for the first moment. Defaults to \code{TRUE}.
#' @param robust.second whether to perform PCA on a robust estimate of the second moment component or not.
#' @param ... additional arguments to pass to the function.
#' @return A list with the following components:
#' \item{\code{A}}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{\code{centroids}}{\code{[K, r]} the centroids of the classes in the projected space.}
#' \item{\code{priors}}{\code{[K]} the prior probabilities of the classes.}
#' \item{\code{ylabs}}{\code{[K]} the unique labels of the classes.}
#' \item{\code{Xr}}{\code{[n, r]} the data projected into the \code{r} dimensions.}
#' \item{\code{cr}}{\code{[K, r]} the centroids projected into the \code{r} dimensions.}
#' \item{\code{second.moment}}{the method used to estimate the second moment.}
#' \item{\code{first.moment}}{the method used to estimate the first moment.}
#' @examples
#' library(spp)
#' data(iris)
#'
#' X <- iris[,1:4]
#' Y <- iris[,5]
#' model <- spp(X, Y, r=3)
#'
#' # use mvi method for projection
#' model <- spp(X, Y, r=2, method='mvi')
#' @export

spp <- function(X, Y, r, method='pca+pls', orth.first=FALSE,
                second.moment.xfm=FALSE, second.moment.xfm.opts=list(),
                first.moment='delta', second.moment='linear', orthogonalize=FALSE,
                robust.first=TRUE, robust.second=FALSE,  ...){
  #class data
  info <- lol.utils.info(X, Y, robust=robust.first)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  if (r > d) {
    stop(sprintf("Embedding dimensions, r=%d, must be lower than the native dimensions, d=%d", r, d))}

  if (first.moment == "delta") {
    first.moment.proj <- lol.utils.deltas(centroids, priors, robust=robust.first)}
  else {first.moment.proj <- array(0, dim=c(d, 0))}
  colnames(first.moment.proj) <- gen.colname('p1',dim(first.moment.proj)[2])

  nv <- r - dim(first.moment.proj)[2]
  if (method %in% c('mvi','xi')){
    A2_all <- project.all(X,Y,r,second.moment.xfm,second.moment.xfm.opts)
    A_all <- cbind(first.moment.proj, A2_all)
    A <- proj.screen(X, Y, A_all, r = r, screen.by = method, orth.first = orth.first)
  }else if(method=='double.proj'){
    A2_all <- project.all(X,Y,r,second.moment.xfm,second.moment.xfm.opts)
    A_all <- cbind(first.moment.proj,A2_all)
    pca.res <- lol.project.pca(A_all,r) # ! use min(n,d) instead of r!
    A <- A_all %*% pca.res$A
  }else if(method == 'pca+pls'){
    lrlda <- lol.project.lrlda(X, Y, r=ceiling(nv/2), xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts, robust=robust.second)
    second.moment.proj_1 <- lrlda$A
    pls.res <- lol.project.pls(X, Y, r=nv-ceiling(nv/2))
    second.moment.proj_2 <- pls.res$A
    # second.moment.proj = array(0, dim=c(d, nv))
    # second.moment.proj[,1] = second.moment.proj_1[,1]
    # for (i in 2:nv) {
    #   if(i%%2 == 0 ){
    #     second.moment.proj[,i] = second.moment.proj_2[,i/2]
    #   }else{
    #     second.moment.proj[,i] = second.moment.proj_1[,(i+1)/2]
    #   }
    # }
    A_all=cbind(second.moment.proj_1, second.moment.proj_2)
    A = cbind(first.moment.proj, A_all)
  }else{
    stop('Unkown method!')
  }
  if (orthogonalize) { A <- qr.Q(qr(A))}
  centroids <- t(centroids)
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs, Xr=lol.embed(X, A),
              cr=lol.embed(centroids, A), second.moment=second.moment, first.moment=first.moment))
}

first.moment <- function(X,Y,r,robust,first.moment = "delta") {
  info <- lol.utils.info(X, Y,r, robust)
  priors <- info$priors
  centroids <- info$centroids
  K <- info$K
  ylabs <- info$ylabs
  n <- info$n
  d <- info$d
  if (r > d) {stop(sprintf("Embedding dimensions, r=%d, must be lower than the native dimensions, d=%d", r, d))}
  if (first.moment == "delta") {
    proj1 <- lol.utils.deltas(centroids, priors)}
  else { proj1 <- array(0, dim=c(d, 0))}
  return(proj1)
}

project.all <- function(X, Y,r,second.moment.xfm=FALSE, second.moment.xfm.opts=list(),...) {
  n <- dim(X)[1]
  d <- dim(X)[2]

  # proj 1: mean difference. #! do not screen!
  # proj1   <- first.moment(X,Y,r,robust=FALSE)
  # colnames(proj1) <- gen.colname('p1',dim(proj1)[2])
  # proj1.r <- first.moment(X,Y,r,robust=TRUE)
  # colnames(proj1.r) <- gen.colname('p1.r',dim(proj1.r)[2])

  # second moment: linear
  lda.res <- lol.project.lrlda(X, Y, r, xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts, robust=FALSE)
  proj2.lda <- lda.res$A
  colnames(proj2.lda) <- gen.colname('p2.lda',r)
  # lrlda.robust <- lol.project.lrlda(X, Y, r, xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts, robust=TRUE)
  # proj2.lrlda.robust <- lrlda.robust$A

  # second moment: Quadratic
  #proj2.quad <- second.moment.quadratic(X,Y,r,second.moment.xfm,second.moment.xfm.opts,robust = FALSE) # use nv instead of r?
  # proj2.quad.robust <- second.moment.quadratic(X,Y,r,second.moment.xfm,second.moment.xfm.opts,robust = TRUE)

  # second moment: pls
  pls.res <- lol.project.pls(X, Y, r) #! use r instead of nv
  proj2.pls <- pls.res$A
  colnames(proj2.pls) <- gen.colname('pls',dim(proj2.pls)[2])

  # PCA
  #pca.res <- lol.project.pca(X,min(n,d)) # ! use min(n,d) instead of r!
  # # proj.pca <- pca.res$A[,1:r]
  # proj.pca <- pca.res$A
  # colnames(proj.pca) <- gen.colname('pca',dim(proj.pca)[2])

  # A_all <- cbind(proj1,proj1.robust,proj2.lrlda, proj2.quad, proj2.pls,proj.pca)
  # A_all <- cbind(proj1, proj2.lda, proj2.pls, proj.pca)
  # A_all <- cbind( proj2.lda, proj2.quad, proj2.pls, proj.pca)
  A_all <- cbind(proj2.lda, proj2.pls)
  # print(dim(A_all))
  return(A_all)
}

proj.screen <- function(X,Y,A_all,r,screen.by,orth.first=FALSE) {
  if (orth.first) {
    A_all <- qr.Q(qr(A_all))
  }
  X_proj <- X %*% A_all
  if (screen.by=='mvi') {
    # screen_result <- screenIID(X_proj,Y,method = 'MV-SIS')
    screen_result <- VariableScreening::screenIID(X_proj,Y,method = 'MV-SIS')
    # cat(screen_result$rank[1:10])
    idx <- screen_result$rank[1:r]
  }else if(screen.by=='xi'){
    n_projs <- dim(X_proj)[2]
    xi_vec <- rep(0, n_projs)
    Y_ <- as.numeric(Y)
    for (j in 1:n_projs) {
      xi_vec[j] <- XICOR::xicor(X_proj[,j],Y_)}
    idx <- order(xi_vec,decreasing = TRUE)[1:r]
  }else if(screen.by=='None'){
    idx <- 1:dim(A_all)[2]
  } else{
    stop(print('Wrong screening method!'))
  }
  # cat(screen.by,':',idx,'\n')
  A <- A_all[,idx]
  return(A)
}


lol.utils.info <- function(X, Y, robust=FALSE, ...) {
  require(robustbase)
  ylabs <- as.vector(sort(unique(Y)))
  dimx <- dim(X)
  n <- dimx[1]; d <- dimx[2]
  K <- length(ylabs)
  # compute the fraction of each class for prior p
  priors = sapply(ylabs, function(y) sum(Y == y)/n)
  if (!robust) {
    centroids <- as.matrix(array(sapply(ylabs, function(y) colMeans(X[Y==y,,drop=FALSE])), dim=c(d, K)))
  } else {
    # robust estimator of the mean is the median
    centroids <- as.matrix(array(sapply(ylabs, function(y) colMedians(X[Y==y,,drop=FALSE])), dim=c(d, K)))
  }
  return(list(n=n, d=d, ylabs=ylabs, priors=priors, centroids=centroids, K=K))
}

