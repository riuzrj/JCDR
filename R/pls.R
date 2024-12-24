#' Partial least-squares based on lolR (Eric)
#' @importFrom pls plsr


spp.project.pls <- function(X, Y, r, ...) {
  info <- spp.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d

  Yh <- spp.utils.ohe(Y)
  tmpData <- data.frame(n=paste("row", 1:nrow(Yh), sep=""))
  tmpData$Y <- Yh
  tmpData$X <- X

  A <- plsr(Y ~ X, data=tmpData, ncomp=r)$projection
  centroids <- t(centroids)
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=spp.embed(X, A), cr=spp.embed(centroids, A)))
}
