#' Residualization function
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param a nxq confounding matrix
#' @param parallel flag for whether to use apply or parApply 
#' @param clust cluster of cores for parallization
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples get.residual(X, Y, a, parallel=F, clust=NULL)

get.residual <- function(X, Y, a, parallel=F, clust=NULL){
  X = as.matrix(X) # make sure it's a matrix
  a = as.matrix(a)
  if(parallel){
    return(get.resid.parallel(X=X,Y=Y,a=a,clust=clust))
  }
  else{
    return(get.resid.apply(X=X,Y=Y,a=a))
  }
}