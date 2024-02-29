#' SIMPLS deflation step
#' @param X nxp data matrix
#' @param V matrix of principal component vectors 
#' @param n number of observations 
#' @param p number of predictors 
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples simpls.deflation.regular(X,V,n,p)

simpls.deflation.regular <- function(X,V,n,p){
  Z <- X%*%V
  P <- crossprod(X,Z)%*%solve(crossprod(Z,Z))
  Q <- diag(x=1,nrow=p) - P%*%tcrossprod(solve(crossprod(P,P)),P)
  return(X%*%Q)
}