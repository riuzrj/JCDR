#' class-conditional PCA based on lolR (Eric)



spp.project.lrlda <- function(X, Y, r, xfm=FALSE, xfm.opts=list(), robust=FALSE, ...) {
  # class data
  classdat <- spp.utils.info(X, Y, robust=robust)
  priors <- classdat$priors; centroids <- t(classdat$centroids)
  K <- classdat$K; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }

  # subtract column means per-class
  Yidx <- sapply(Y, function(y) which(ylabs == y))
  # form class-conditional data matrix
  Xt <- X - centroids[Yidx,]
  # compute the standard projection but with the pre-centered data.
  X.decomp <- spp.utils.decomp(Xt, xfm=xfm, xfm.opts=xfm.opts, ncomp=r, robust=robust)

  return(list(A=X.decomp$comp, d=X.decomp$val, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=spp.embed(X, X.decomp$comp), cr=spp.embed(centroids, X.decomp$comp)))
}



spp.utils.decomp <- function(X, xfm=FALSE, xfm.opts=list(), ncomp=0, t=.05, robust=FALSE) {
  n <- nrow(X)
  d <- ncol(X)
  # scale if desired before taking SVD
  for (i in 1:length(xfm)) {
    sc <- xfm[i]
    if (!(i %in% names(xfm.opts))) {
      xfm.opts[[i]] <- list()
    }
    if (sc == 'unit') {
      X <- do.call(scale, c(list(X), xfm.opts[[i]]))
    } else if (sc == 'log') {
      X <- do.call(log, c(list(X), xfm.opts[[i]]))
    } else if (sc == 'rank') {
      X <- apply(X, c(2), function(x) {do.call(rank, c(list(x), xfm.opts[[i]]))})
    }
  }
  # take svd
  decomp <- list()
  if (robust) {
    eigenX <- eigen(covRob(X, estim='weighted')$cov)
    decomp$comp <- eigenX$vectors[, 1:ncomp, drop=FALSE]
    decomp$val <- eigenX$values[1:ncomp]
  } else if (ncomp > t*d | ncomp >= d) {
    svdX <- svd(X, nu=0, nv=ncomp)
    decomp$comp <- svdX$v
    decomp$val <- svdX$d
  } else {
    svdX <- irlba(X, nu=0, nv=ncomp)
    decomp$comp <- svdX$v
    decomp$val <- svdX$d
  }
  return(decomp)
}
