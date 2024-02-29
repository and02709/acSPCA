#' Eigendecomposition using residualization process
#' @param X.star Data matrix residualized with respect to confounding variables
#' @param X.tilde Data matrix residualized with respect to response variable
#' @param A kernel matrix of confounders
#' @param L kernel matrix of response
#' @param H centering matrix 
#' @param lambda penalty term for confounding variable adjustment 
#' @param npc number of principl components
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples reg.mat.eig.sym.decomp(X.star,X.tilde,A,L,H,lambda,npc)

reg.mat.eig.sym.decomp <- function(X.star,X.tilde,A,L,H,lambda,npc){
  return(RSpectra::eigs_sym(crossprod(X.star,H)%*%L%*%H%*%X.star - lambda*crossprod(X.tilde,H)%*%A%*%H%*%X.tilde, k=npc, which="LA"))
}