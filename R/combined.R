#‘ A function that performs the combination
#‘
#‘
#' @export

project.all <- function(X, Y,r,second.moment.xfm=FALSE, second.moment.xfm.opts=list(),...) {
  n <- dim(X)[1]
  d <- dim(X)[2]

  # second moment: linear
  lda.res <- spp.project.lrlda(X, Y, r, xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts, robust=FALSE)
  proj2.lda <- lda.res$A
  colnames(proj2.lda) <- gen.colname('p2.lda',r)

  # second moment: pls
  pls.res <- spp.project.pls(X, Y, r) #! use r instead of nv
  proj2.pls <- pls.res$A
  colnames(proj2.pls) <- gen.colname('pls',dim(proj2.pls)[2])

  A_all <- cbind(proj2.lda, proj2.pls)
  # print(dim(A_all))
  return(A_all)
}
