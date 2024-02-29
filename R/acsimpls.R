#' SIMPLS algorithm included adjustment for confounding variables
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param A nxq confounding variable matrix
#' @param resp.kernel kernel of response vector which can be linear or delta 
#' @param conf.kernel kernel of confounding matrix which can be linear or gaussian 
#' @param bandwidth gaussian kernel argument 
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples acsimpls(X,Y,A,lambda,npc,resp.kernel="linear", conf.kernel="linear", bandwidth=NULL)


acsimpls <- function(X,Y,A,lambda,npc,resp.kernel=c("linear","delta"), conf.kernel=c("linear", "gaussian"), bandwidth=NULL){
  nobs <- nrow(X)
  p <- ncol(X)
  Ka <- calkernel(A, conf.kernel, bandwidth)
  Ky <- respkernel(Y, dim(Y)[[1]], resp.kernel)
  H <- diag(x=1, nrow=nobs) - 1/nobs*rep(1,nobs)%*%t(rep(1,nobs))
  
  V <- matrix(0, nrow=p, ncol=npc)
  Xi <- X
  
  for(i in 1:npc){
    spca.obj <- single.vector.decomp(X=Xi,H=H,Ky=Ky,Ka=Ka,lambda=lambda)
    V[,i] <- spca.obj$vectors
    if(i==npc) break
    Xi <- deflate(X=X,V=V[,1:i],p=p)
  }
  return(V)
}