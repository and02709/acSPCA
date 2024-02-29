#' SIMPLS deflation step using Rfast
#' @param X nxp data matrix
#' @param V matrix of principal component vectors 
#' @param n number of observations 
#' @param p number of predictors 
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples simpls.deflation.rfast(X,V,n,p)

simpls.deflation.rfast <- function(X,V,n,p){
  Z <- X%*%V
  #P <- mat.mult(Crossprod(X,Z),spdinv(Crossprod(Z,Z)))
  P <- Crossprod(X,Z)%*%spdinv(crossprod(Z,Z))
  #Q <- Diag.matrix(p,1) - mat.mult(P,Tcrossprod(spdinv(Crossprod(P,P)),P))
  Q <- diag(x=1,nrow=p) - P%*%tcrossprod(spdinv(crossprod(P,P)),P)
  return(mat.mult(X,Q))
}