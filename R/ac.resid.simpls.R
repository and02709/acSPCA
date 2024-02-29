#' Residualized adjusted for confounding SIMPLS
#' @param X nxp data matrix
#' @param Y nx1 response matrix
#' @param a nxq confounding variable matrix
#' @param A kernel matrix of confounders
#' @param L kernel matrix of response
#' @param H centering matrix 
#' @param lambda penalty term for confounding variable adjustment 
#' @param npc number of principl components
#' @param parallel parallelization of residualization 
#' @param clust cluster of cores for parallel procedure
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples ac.resid.simpls(X, Y, a, L, A, H, lambda, npc, parallel=F, clust=NULL)

ac.resid.simpls <- function(X, Y, a, L, A, H, lambda, npc, parallel=F, clust=NULL){
  n <- dim(X)[[1]]
  p <- dim(X)[[2]]
  V <- matrix(0, nrow = p, ncol=npc)
  Xi <- X
  for(i in 1:npc){
    X.temp <- get.residual(X=Xi, Y=Y, a=a, parallel=parallel, clust=clust)
    X.star <- matrix(X.temp$X.star, nrow = n, ncol=p)
    X.tilde <- matrix(X.temp$X.tilde, nrow = n, ncol=p)
    V[,i] <- reg.mat.eig.sym.decomp(X.star=X.star,X.tilde=X.tilde,A=A,L=L,H=H,lambda=lambda,npc=1)$vectors
    if(i==npc) break
    Xi <- simpls.deflation.regular(X=X,V=V[,1:i],n=n,p=p)
  }
  return(V)
}