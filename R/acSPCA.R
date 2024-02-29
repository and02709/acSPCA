#' Supervised Adjusted for Confounding Principal Component Analysis
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param A nxq confounding variable matrix
#' @param lambda penalty to adjust for confounding
#' @param npc number of desired principal component
#' @param resp.kernel kernel of response vector which can be linear or delta 
#' @param conf.kernel kernel of confounding matrix which can be linear or gaussian 
#' @param bandwidth gaussian kernel argument 
#' @keywords Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples acSSPCA(X, Y, A, X4A=NULL, X4Y=NULL, c1=NULL, c2, v_ini, resp.kernel=c("linear","delta"), conf.kernel=c("linear", "gaussian"), bandwidth=NULL, maxiter=50, delta=10^-8, filter=T, minmaxsep=1e-3)

acSPCA <- function(X,Y,A,lambda,npc,resp.kernel=c("linear","delta"), conf.kernel=c("linear", "gaussian"), bandwidth=NULL){
  
  nobs <- dim(X)[[1]]
  
  #missing data
  X[is.na(X)] <- mean(X, na.rm=T)
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  A[is.na(A)] <- mean(A, na.rm=T)
  
  Ka <- calkernel(A, conf.kernel, bandwidth)
  Ky <- respkernel(Y, dim(Y)[[1]], resp.kernel)
  H <- diag(x=1, nrow=nobs) - 1/nobs*rep(1,nobs)%*%t(rep(1,nobs))
  
  return(eigs_sym(mat.mult(mat.mult(mat.mult(Crossprod(X,H),Ky),H),X) - lambda*mat.mult(mat.mult(mat.mult(Crossprod(X,H),Ka),H),X), k=npc, which="LA"))
}