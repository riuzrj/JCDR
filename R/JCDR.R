#' Joint Combined Dimensionality Reduction (JCDR)
#'
#' A function that performs supervised PCA-PLS.
#'
#'
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
#' library(JCDR)
#' data(iris)
#'
#' X <- as.matrix(iris[,1:4])
#' Y <- iris[,5]
#' model <- spp(X, Y, r=4)
#'
#' # use mvi method for projection
#' model <- spp(X, Y, r=4, method='mvi')
#' @export

spp <- function(X, Y, r, method='pca+pls', orth.first=FALSE,
                second.moment.xfm=FALSE, second.moment.xfm.opts=list(),
                first.moment='delta', second.moment='linear', orthogonalize=FALSE,
                robust.first=TRUE, robust.second=FALSE,  ...){
  #class data
  info <- spp.utils.info(X, Y, robust=robust.first)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  if (r > d) {
    stop(sprintf("Embedding dimensions, r=%d, must be lower than the native dimensions, d=%d", r, d))}

  if (first.moment == "delta") {
    first.moment.proj <- spp.utils.deltas(centroids, priors, robust=robust.first)}
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
    pca.res <- spp.project.lrlda(A_all,r) # ! use min(n,d) instead of r!
    A <- A_all %*% pca.res$A
  }else if(method == 'pca+pls'){
    lrlda <- spp.project.lrlda(X, Y, r=ceiling(nv/2), xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts, robust=robust.second)
    second.moment.proj_1 <- lrlda$A
    pls.res <- spp.project.pls(X, Y, r=nv-ceiling(nv/2))
    second.moment.proj_2 <- pls.res$A
    A_all=cbind(second.moment.proj_1, second.moment.proj_2)
    A = cbind(first.moment.proj, A_all)
  }else{
    stop('Unkown method!')
  }
  if (orthogonalize) { A <- qr.Q(qr(A))}
  centroids <- t(centroids)
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs, Xr=spp.embed(X, A),
              cr=spp.embed(centroids, A), second.moment=second.moment, first.moment=first.moment))
}



gen.colname <- function(model_name,r){
  name_vec <- c()
  for (i in 1:r) {
    name_vec <- c(name_vec,paste(model_name,toString(i),sep='.'))  }
  return(name_vec)
}

spp.embed <- function(X, A, ...) {
  return(X %*% A)
}

# first.moment <- function(X,Y,r,robust,first.moment = "delta") {
#   info <- spp.utils.info(X, Y,r, robust)
#   priors <- info$priors
#   centroids <- info$centroids
#   K <- info$K
#   ylabs <- info$ylabs
#   n <- info$n
#   d <- info$d
#   if (r > d) {stop(sprintf("Embedding dimensions, r=%d, must be lower than the native dimensions, d=%d", r, d))}
#   if (first.moment == "delta") {
#     proj1 <- spp.utils.deltas(centroids, priors)}
#   else { proj1 <- array(0, dim=c(d, 0))}
#   return(proj1)
# }








