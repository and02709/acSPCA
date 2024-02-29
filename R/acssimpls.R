#' Sparse SIMPLS algorithm using adjusted for confounding variables
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param A nxq confounding variable matrix
#' @param lambda penalty for adjusting for confounding 
#' @param npc number of principal components
#' @param resp.kernel kernel of response vector which can be linear or delta 
#' @param conf.kernel kernel of confounding matrix which can be linear or gaussian 
#' @param c1 norm of quadratic predictor and confounding matrix product
#' @param c2 sum of absolute value loadings for desired sparse vector 
#' @param bandwidth gaussian kernel argument 
#' @param maxiter maximum number of iterations for the algorithm 
#' @param delta offset parameter
#' @param filter parameter to determine filtering
#' @param minmaxsep determines the acceptable algorithm termination acceptable different
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples acssimpls(X,Y,A,lambda,npc,resp.kernel="linear", conf.kernel="linear", bandwidth=NULL, c1=NULL, c2=NULL, maxiter=50, delta=10^-8, filter=T, minmaxsep=1e-3)

acssimpls <- function(X,Y,A,lambda,npc,resp.kernel=c("linear","delta"), conf.kernel=c("linear", "gaussian"), bandwidth=NULL, 
                      c1=NULL, c2=NULL, maxiter=50, delta=10^-8, filter=T, minmaxsep=1e-3){
  nobs <- nrow(X)
  p <- ncol(X)
  Ka <- calkernel(A, conf.kernel, bandwidth)
  Ky <- respkernel(Y, dim(Y)[[1]], resp.kernel)
  H <- diag(x=1, nrow=nobs) - 1/nobs*rep(1,nobs)%*%t(rep(1,nobs))
  
  V <- matrix(0, nrow=p, ncol=npc)
  Xi <- X
  
  for(i in 1:npc){
    spca.obj <- single.vector.decomp(X=Xi,H=H,Ky=Ky,Ka=Ka,lambda=lambda)
    vtemp <- spca.obj$vectors
    sspca.obj <- acSSPCA(X=X, Y=Y, A=A, X4A=NULL, X4Y=NULL, c1=NULL, c2=c2, v_ini=vtemp, resp.kernel=resp.kernel, conf.kernel=conf.kernel, 
                         bandwidth=bandwidth, maxiter=maxiter, delta=delta, filter=filter, minmaxsep=minmaxsep)
    vs <- sspca.obj$v
    V[,i] <- vs
    if(i==npc) break
    Xi <- deflate(X=X,V=V[,1:i],p=p)
  }
  return(V)
}
