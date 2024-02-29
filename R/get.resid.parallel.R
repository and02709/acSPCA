#' Uses parApply function to perform residual calculation
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param a nxq confounding matrix
#' @param clust cluster of cores for parallelization
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples get.resid.parallel(X, Y, a, clust)

get.resid.parallel <- function(X, Y, a, clust){
  temp <- parApply(cl=clust,X=X,2,res.calc,Y=Y,a=a)
  X.star <- unlist(lapply(X=temp, function(x) return(x$star)))
  X.tilde <- unlist(lapply(X=temp, function(x) return(x$tilde)))
  return(list(X.tilde=X.tilde,X.star=X.star))
}