#' Deflation of the X matrix in the SIMPLS algorithm
#' @param X nxp predictor matrix
#' @param V matrix containing principal component vectors
#' @param p number of predictors in data matrix X
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples deflate(X,V,p)

deflate <- function(X,V,p){
  # Encode vectors to data matrix
  Z <- X%*%as.matrix(V)
  
  # Generate spanning vector
  #P <- t(X)%*%Z%*%solve(t(Z)%*%Z)
  P <- crossprod(X,Z%*%solve(crossprod(Z,Z)))
  
  # Generate orthogonal projection matrix
  #Q <- diag(x=1,nrow = p) - P%*%solve(t(P)%*%P)%*%t(P)
  Q <- diag(x=1,nrow = p) - P%*%tcrossprod(solve(crossprod(P,P)),P)
  
  return(X%*%Q)
}