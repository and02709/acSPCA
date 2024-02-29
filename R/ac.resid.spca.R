#' Residualized adjusted for confounding supervised principal component analysis
#' @param X nxp data matrix
#' @param Y nx1 response matrix
#' @param a nxq confounding variable matrix
#' @param A kernel matrix of confounders
#' @param L kernel matrix of response
#' @param H centering matrix 
#' @param lambda penalty term for confounding variable adjustment 
#' @param npc number of principl components
#' @param fast.mat flag for whether to use Rfast coding
#' @param parallel parallelization of residualization 
#' @param clust cluster of cores for parallel procedure
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples ac.resid.spca(X, Y, a, L, A, H, lambda, npc, fast.mat=F, parallel=F, clust=NULL)

ac.resid.spca <- function(X, Y, a, L, A, H, lambda, npc, fast.mat=F, parallel=F, clust=NULL){
  # Residualize the data matrix X
  n <- nrow(X)
  p <- ncol(X)
  resid.dat = get.residual(X=X, Y=Y, a=a, parallel=parallel,clust=clust)
  X.star = matrix(resid.dat$X.star,nrow = n,ncol = p)
  X.tilde = matrix(resid.dat$X.tilde,nrow = n,ncol = p)
  return(reg.mat.eig.sym.decomp(X.star=X.star,X.tilde=X.tilde,A=A,
                                  L=L,H=H,lambda=lambda,npc=npc)$vectors)

}