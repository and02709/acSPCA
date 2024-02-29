#' Calculates the Psi matix Psi=Delta^t HX.
#' @param X nxp data matrix.  Typically should be centered
#' @param L nxn response kernel matrix
#' @param n number of responses in the data
#' @param p number of predictors in the data
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples get.psi(X,L,n,p)

get.psi <- function(X,L,n,p){
  # Centering matrix
  H <- diag(x=1, nrow=n) - 1/n*rep(1,n)%*%t(rep(1,n))
  # Eigendecomposition of L
  Eigendecomp <- eigen(L)
  # Need to retain eigenvectors to reconstruct Delta
  U <- Eigendecomp$vectors
  # Eigenvalues of kernel response matrix L
  EV <- Eigendecomp$values
  # Generation of diagonal matrix of square root eigenvalues of L
  Sigmat <- diag(sqrt(zapsmall(EV)))
  # Generation of Delta such that t(Delta)%*%Delta=L
  #Delta <- Rfast::Tcrossprod(Rfast::mat.mult(U, Sigmat), U)
  Delta <- tcrossprod(U%*%Sigmat,U)
  # Generation of Psi matrix
  #return(Rfast::Crossprod(Delta, Rfast::mat.mult(H,X)))
  return(crossprod(Delta,H%*%X))
}