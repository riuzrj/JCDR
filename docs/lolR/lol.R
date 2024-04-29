#' Linear Optimal Low-Rank Projection (LOL)
#'
#' A function for implementing the Linear Optimal Low-Rank Projection (LOL) Algorithm. This algorithm allows users to find an optimal
#' projection from `d` to `r` dimensions, where `r << d`, by combining information from the first and second moments in thet data.
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
#' \itemize{
#' \item{\code{'linear'} performs PCA on the class-conditional data to capture the second moment, retaining the vectors with the top singular values.   Transform options for \code{second.moment.xfm} and arguments in \code{second.moment.opts} should be in accordance with the trailing arguments for \link{lol.project.lrlda}.}
#' \item{\code{'quadratic'} performs PCA on the data for each class separately to capture the second moment, retaining the vectors with the top singular values from each class's PCA. Transform options for \code{second.moment.xfm} and arguments in \code{second.moment.opts} should be in accordance with the trailing arguments for \link{lol.project.pca}.}
#' \item{\code{'pls'} performs PLS on the data to capture the second moment, retaining the vectors that maximize the correlation between the different classes. Transform options for \code{second.moment.xfm} and arguments in \code{second.moment.opts} should be in accordance with the trailing arguments for \link{lol.project.pls}.}
#' \item{\code{FALSE} do not capture the second moment.}
#' }
#' @param orthogonalize whether to orthogonalize the projection matrix. Defaults to \code{FALSE}.
#' @param second.moment.xfm whether to use extraneous options in estimation of the second moment component. The transforms specified should be a numbered list of transforms you wish to apply, and will be applied in accordance with \code{second.moment}.
#' @param second.moment.xfm.opts optional arguments to pass to the \code{second.moment.xfm} option specified. Should be a numbered list of lists, where \code{second.moment.xfm.opts[[i]]} corresponds to the optional arguments for \code{second.moment.xfm[[i]]}.
#' Defaults to the default options for each transform scheme.
#' @param robust.first whether to perform PCA on a robust estimate of the first moment component or not. A robust estimate corresponds to usage of medians. Defaults to \code{TRUE}.
#' @param robust.second whether to perform PCA on a robust estimate of the second moment component or not. A robust estimate corresponds to usage of a robust covariance matrix, which requires \code{d < n}. Defaults to \code{FALSE}.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{A}}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{\code{ylabs}}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{\code{centroids}}{\code{[K, d]} centroid matrix of the \code{K} unique, ordered classes in native \code{d} dimensions.}
#' \item{\code{priors}}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{\code{Xr}}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' \item{\code{cr}}{\code{[K, r]} the \code{K} centroids in reduced dimensionality \code{r}.}
#' \item{\code{second.moment}}{the method used to estimate the second moment.}
#' \item{\code{first.moment}}{the method used to estimate the first moment.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("lol", package = "lolR")}
#'
#' @author Eric Bridgeford

#' @references Joshua T. Vogelstein, et al. "Supervised Dimensionality Reduction for Big Data" arXiv (2020).
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.lol(X=X, Y=Y, r=5)  # use lol to project into 5 dimensions
#'
#' # use lol to project into 5 dimensions, and produce an orthogonal basis for the projection matrix
#' model <- lol.project.lol(X=X, Y=Y, r=5, orthogonalize=TRUE)
#'
#' # use LRQDA to estimate the second moment by performing PCA on each class
#' model <- lol.project.lol(X=X, Y=Y, r=5, second.moment='quadratic')
#'
#' # use PLS to estimate the second moment
#' model <- lol.project.lol(X=X, Y=Y, r=5, second.moment='pls')
#'
#' # use LRLDA to estimate the second moment, and apply a unit transformation
#' # (according to scale function) with no centering
#' model <- lol.project.lol(X=X, Y=Y, r=5, second.moment='linear', second.moment.xfm='unit',
#'                          second.moment.xfm.opts=list(center=FALSE))
#' @export
lol.project.lol <- function(X, Y, r, method='lol',fix.proj1=FALSE, orth.first=FALSE,second.moment.xfm=FALSE, second.moment.xfm.opts=list(),
                            first.moment='delta', second.moment='linear', orthogonalize=FALSE,
                            robust.first=TRUE, robust.second=FALSE,  ...) {
    
    # class data
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
    if (second.moment == "linear" & nv > 0) {
        lrlda <- lol.project.lrlda(X, Y, r=nv, xfm=second.moment.xfm,xfm.opts=second.moment.xfm.opts, robust=robust.second)
        second.moment.proj <- lrlda$A}
    else if (second.moment == "quadratic" & nv > 0) {
        Aclass <- array(0, dim=c(d, 0))  # the class-wise egvecs
        vclass <- c()  # the class-wise egvals
        for (ylab in ylabs) {
            Xclass = X[Y == ylab,]
            obj <- lol.project.pca(Xclass, r=nv, xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts, robust=robust.second)
            Aclass <- cbind(Aclass, obj$A)
            vclass <- c(vclass, obj$d[1:nv])}
        # take the nv from the A computed for each class using the nv with the top eigenvalues from Aclass
        second.moment.proj <- Aclass[, sort(vclass, index.return=TRUE, decreasing=TRUE)$ix[1:nv]]}
    else if (second.moment == "pls" & nv > 0) {
        pls.res <- lol.project.pls(X, Y, nv)
        second.moment.proj <- pls.res$A}
    else {second.moment.proj <- array(0, dim=c(d, 0))}
    
    if (method=='lol') {
        A <- cbind(first.moment.proj, second.moment.proj)}
    else if(method!='lol' & fix.proj1==TRUE){
        if (nv>0) {
            A2_all <- project.all(X,Y,r,second.moment.xfm,second.moment.xfm.opts)
            # keep second.moment.projs by number of r=nv .
            A2 <- proj.screen(X,Y,A2_all,r = nv, screen.by = screen.by, orth.first = orth.first) }
        else{ A2 <- array(0, dim=c(d, 0))}
        A <- cbind(first.moment.proj,A2) # keep first.moment.proj
    }else if(method %in% c('mvi','xi')){
        #print('mvi')
        A2_all <- project.all(X,Y,r,second.moment.xfm,second.moment.xfm.opts)
        A_all <- cbind(first.moment.proj, A2_all)
        A <- proj.screen(X, Y, A_all, r = r, screen.by = method, orth.first = orth.first)
        #print('xi')
    }else if(method=='double.proj'){
        A2_all <- project.all(X,Y,r,second.moment.xfm,second.moment.xfm.opts)
        A_all <- cbind(first.moment.proj,A2_all)       
        pca.res <- lol.project.pca(A_all,r) # ! use min(n,d) instead of r!
        A <- A_all %*% pca.res$A
    }else if(method=='lol1'){
        A_all <- project.some1(X, Y,r,second.moment.xfm, second.moment.xfm.opts)
        A <- cbind(first.moment.proj, A_all)
    }else if(method=='QOQ1'){
        A_all <- project.some2(X, Y,r,second.moment.xfm, second.moment.xfm.opts)
        A <- cbind(first.moment.proj, A_all)
    }else if(method == 'pca+pls'){
      #print(nv)
      lrlda <- lol.project.lrlda(X, Y, r=ceiling(nv/2), xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts, robust=robust.second)
      
      second.moment.proj_1 <- lrlda$A
      #print(dim(second.moment.proj_1))
      #print(second.moment.proj_1)
      pls.res <- lol.project.pls(X, Y, r=nv-ceiling(nv/2))
      second.moment.proj_2 <- pls.res$A
      
      #print(dim(second.moment.proj_2))
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
      #print(A_all)
      A = cbind(first.moment.proj, A_all)
        #A_all <- project.some3(X, Y,nv,second.moment.xfm, second.moment.xfm.opts,robust=robust)
        
        #print(A_all)
        #A = cbind(first.moment.proj, A_all)
        
    }else if(method == 'select_pca+pls'){
      
      lrlda <- lol.project.lrlda(X, Y, r=nv, xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts, robust=robust.second)
      second.moment.proj_1 <- lrlda$A
      
      pls.res <- lol.project.pls(X, Y, r=nv)
      second.moment.proj_2 <- pls.res$A
      
      second.moment.proj = cbind(second.moment.proj_1, second.moment.proj_2)
      
      A = cbind(first.moment, second.moment.proj)
      
      
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


second.moment.quadratic <- function(X,Y,r,second.moment.xfm,second.moment.xfm.opts,robust) {
    info <- lol.utils.info(X, Y, robust)
    priors <- info$priors
    centroids <- info$centroids
    K <- info$K
    ylabs <- info$ylabs
    n <- info$n
    d <- info$d
    
    Aclass <- array(0, dim=c(d, 0))  # the class-wise egvecs
    vclass <- c()  # the class-wise egvals
    for (ylab in ylabs) {
        Xclass = X[Y == ylab,]
        obj <- lol.project.pca(Xclass, r=r,xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts,robust=robust)
        Aclass <- cbind(Aclass, obj$A)
        vclass <- c(vclass, obj$d[1:r])}
    proj2.quadratic <- Aclass[, sort(vclass, index.return=TRUE, decreasing=TRUE)$ix[1:r]]
    return(proj2.quadratic)
}

gen.colname <- function(model_name,r){
    name_vec <- c()
    for (i in 1:r) {
        name_vec <- c(name_vec,paste(model_name,toString(i),sep='.'))  }
    return(name_vec)
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
    
# project.some1 <- function(X, Y,r,second.moment.xfm=FALSE, second.moment.xfm.opts=list(),..){
#     n <- dim(X)[1]
#     d <- dim(X)[2]
#     # linear + PLS
#     lda.res <- lol.project.lrlda(X, Y, r, xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts, robust=FALSE)
#     proj2.lda <- lda.res$A
#     colnames(proj2.lda) <- gen.colname('p2.lda',r)
#     
#     pls.res <- lol.project.pls(X, Y, r) #! use r instead of nv
#     proj2.pls <- pls.res$A
#     colnames(proj2.pls) <- gen.colname('pls',dim(proj2.pls)[2])
#     
#     A_all <- cbind( proj2.lda, proj2.pls)
#     # print(dim(A_all))
#     return(A_all)
# }
# 
# project.some2 <- function(X, Y,r,second.moment.xfm=FALSE, second.moment.xfm.opts=list(),..){
#     n <- dim(X)[1]
#     d <- dim(X)[2]
#     # Quadratic + PLS
#     proj2.quad <- second.moment.quadratic(X,Y,r,second.moment.xfm,second.moment.xfm.opts,robust = FALSE)
#     
#     pls.res <- lol.project.pls(X, Y, r) #! use r instead of nv
#     proj2.pls <- pls.res$A
#     colnames(proj2.pls) <- gen.colname('pls',dim(proj2.pls)[2])
#     
#     A_all <- cbind( proj2.quad, proj2.pls)
#     # print(dim(A_all))
#     return(A_all)
# }

# project.some3 <- function(X, Y,nv,second.moment.xfm=FALSE, second.moment.xfm.opts=list(),robust=robust){
#   n <- dim(X)[1]
#   d <- dim(X)[2]
#   # PCA + PLS
#   print('before')
#   lrlda <- lol.project.lrlda(X, Y, r=ceiling(nv/2), xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts, robust=robust.second)
#   print('after')
#   second.moment.proj_1 <- lrlda$A
#   pls.res <- lol.project.pls(X, Y, r=nv-ceiling(nv/2))
#   second.moment.proj_2 <- pls.res$A
#   A_all = array(0, dim=c(d, nv))
#   for (i in 1:nv) {
#     if(i%%2 == 0 ){
#       A_all[,i] = second.moment.proj_2[,2*i] 
#     }else{
#       A_all[,i] = second.moment.proj_1[,2*i-1]
#     }
#   }
#   print('a')
# 
#   colnames(A_all) <- gen.colname('pca+pls',dim(A_all)[2])
#   
#   
#   # print(dim(A_all))
#   return(A_all)
# }