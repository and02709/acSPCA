#' Deflation of the X matrix in the SIMPLS algorithm
#' @param X nxp predictor matrix
#' @param V matrix containing principal component vectors
#' @param p number of predictors in data matrix X
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples deflate(X,V,p)

deflate <- function(X,V,p){
  # Encode vectors to data matrix
  Z <- Rfast::mat.mult(X,as.matrix(V))
  
  # Generate spanning vector
  #P <- t(X)%*%Z%*%solve(t(Z)%*%Z)
  P <- Rfast::mat.mult(Rfast::Crossprod(X,Z),Rfast::spdinv(Rfast::Crossprod(Z,Z)))
  
  # Generate orthogonal projection matrix
  #Q <- diag(x=1,nrow = p) - P%*%solve(t(P)%*%P)%*%t(P)
  Q <- diag(x=1,nrow = p) - Rfast::mat.mult(P,Rfast::Tcrossprod(Rfast::spdinv(Rfast::Crossprod(P,P)),P))
  
  return(Rfast::mat.mult(X,Q))
}