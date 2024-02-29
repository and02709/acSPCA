#' Orthogonalize design matrix to extracted principal component vector
#' @param X nxp predictor matrix
#' @param v principal component vectors
#' @param nobs number of responses
#' @param p number of predictors
#' @keywords Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples sspca.orthog(X,v,nobs,p)

sspca.orthog <- function(X,v,nobs,p){
  Q <- diag(1,p) - tcrossprod(v,v)
  return(X%*%Q)
}