#' Sparse Supervised Adjusted for Confounding Principal Component Analysis
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param A nxq confounding variable matrix
#' @param lambda penalty for adjusting for confounding 
#' @param npc number of principal components
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
#' @examples acSPCAomni <- function(X,Y,A,lambda,npc,resp.kernel=c("linear","delta"), conf.kernel=c("linear", "gaussian"), bandwidth=NULL, c1=NULL, c2=NULL, maxiter=50, delta=10^-8, filter=T, minmaxsep=1e-3)

acSPCAomni <- function(X,Y,A,lambda,npc,resp.kernel=c("linear","delta"), 
                       conf.kernel=c("linear", "gaussian"), bandwidth=NULL, 
                       c1=NULL, c2=NULL, maxiter=50, delta=10^-8, 
                       filter=T, minmaxsep=1e-3){
  
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
  }
  #if A is a vector, change it to a matrix
  if (is.null(dim(A))){
    A <- matrix(A, ncol=1)
  }
  
  nobs <- dim(X)[[1]]
  
  #check whether a whole row in X is missing
  Xmis <- apply(X, 1, function(row){sum(!is.na(row))})
  if (sum(Xmis==0)){
    stop(paste("The following samples in X is missing, please remove them in X and Y: rows ", paste(which(Xmis==0), collapse =  " "), sep=""))
  }
  
  #check whether the number of samples in X, Y, and A match
  if (dim(X)[1]!=dim(Y)[1]){
    stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
  }
  
  if (dim(X)[1]!=dim(A)[1]){
    stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
  }
  
  if (dim(Y)[1]!=dim(A)[1]){
    stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
  }
  
  #missing data
  X[is.na(X)] <- mean(X, na.rm=T)
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  A[is.na(A)] <- mean(A, na.rm=T)
  
  if(is.null(c2)){
    spca.obj <- acSPCA(X=X,Y=Y,A=A,lambda=lambda,npc=npc,resp.kernel=resp.kernel, conf.kernel=conf.kernel, bandwidth=bandwidth)
    return(list(values=spca.obj$values, vectors=spca.obj$vectors))
  } else{
    sspca.obj <- acSSPCAmulti(X=X,Y=Y,A=A,lambda=lambda,npc=npc,resp.kernel=resp.kernel, conf.kernel=conf.kernel, bandwidth=bandwidth, 
                              c1=c1, c2=c2, maxiter=maxiter, delta=delta, filter=filter, minmaxsep=minmaxsep)
    return(list(values=sspca.obj$values, vectors=sspca.obj$vectors))
  }
}