#' A function that performs the marginal combined dimension reduction (MCDR)
#'
#' @export

proj.screen <- function(X,Y,A_all,r,screen.by,orth.first=FALSE) {
  if (orth.first) {
    A_all <- qr.Q(qr(A_all))
  }
  X_proj <- X %*% A_all
  if (screen.by=='mvi') {
    screen_result <- VariableScreening::screenIID(X_proj,Y,method = 'MV-SIS')
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
