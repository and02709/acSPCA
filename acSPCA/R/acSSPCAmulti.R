#' Sparse Supervised Adjusted for Confounding Principal Component Analysis allowing multiple sparse principal components
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param A nxq confounding variable matrix
#' @param c1 norm of quadratic predictor and confounding matrix product
#' @param c2 sum of absolute value loadings for desired sparse vector 
#' @param resp.kernel kernel of response vector which can be linear or delta 
#' @param conf.kernel kernel of confounding matrix which can be linear or gaussian 
#' @param bandwidth gaussian kernel argument 
#' @param maxiter maximum number of iterations for the algorithm 
#' @param delta offset parameter
#' @param filter parameter to determine filtering
#' @param minmaxsep determines the acceptable algorithm termination acceptable different
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples acSSPCAmulti(X,Y,A,lambda,npc,resp.kernel=c("linear","delta"), conf.kernel=c("linear", "gaussian"), bandwidth=NULL, c1=NULL, c2=NULL, maxiter=50, delta=10^-8, filter=T, minmaxsep=1e-3)


acSSPCAmulti <- function(X,Y,A,lambda,npc,resp.kernel=c("linear","delta"), conf.kernel=c("linear", "gaussian"), bandwidth=NULL, 
                         c1=NULL, c2=NULL, maxiter=50, delta=10^-8, filter=T, minmaxsep=1e-3){
  nobs <- nrow(X)
  p <- ncol(X)
  V <- matrix(0, nrow=p, ncol=npc)
  vals <- rep(0,npc)
  
  #missing data
  X[is.na(X)] <- mean(X, na.rm=T)
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  A[is.na(A)] <- mean(A, na.rm=T)
  
  
  Xi <- X
  for(i in 1:npc){
    spca.obj <- acSPCA(X=Xi,Y=Y,A=A,lambda=lambda,npc=1,resp.kernel=resp.kernel,conf.kernel=conf.kernel,bandwidth = bandwidth)
    vals[i] <- spca.obj$values
    vtemp <- spca.obj$vectors
    sspca.obj <- acSSPCA(X=X, Y=Y, A=A, X4A=NULL, X4Y=NULL, c1=NULL, c2=c2, 
                         v_ini=vtemp, resp.kernel=resp.kernel, 
                         conf.kernel=conf.kernel, 
                         bandwidth=bandwidth, maxiter=maxiter, 
                         delta=delta, filter=filter, minmaxsep=minmaxsep)
    vs <- sspca.obj$v
    V[,i] <- vs
    Xi <- sspca.orthog(X=Xi,v=vs,nobs=nobs,p=p)
  }
  return(list(values=vals,vectors=V))
}