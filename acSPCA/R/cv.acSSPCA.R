#' CV function for determining sparsity penalty c2
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param A nxq confounding variable matrix
#' @param lambda penalty for adjusting for confounding 
#' @param npc number of principal components
#' @param c1 norm of quadratic predictor and confounding matrix product
#' @param c2 vector of values given as sum of absolute value loadings for desired sparse vector 
#' @param resp.kernel kernel of response vector which can be linear or delta 
#' @param conf.kernel kernel of confounding matrix which can be linear or gaussian 
#' @param bandwidth gaussian kernel argument 
#' @param maxiter maximum number of iterations for the algorithm 
#' @param delta offset parameter
#' @param filter parameter to determine filtering
#' @param minmaxsep determines the acceptable algorithm termination acceptable different 
#' @param n.cores number of cores for parallel option 
#' @param parallel flag for whether or not parallelization should be used
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples cv.acSSPCA <- function(X,Y,A,lambda,npc,resp.kernel=c("linear","delta"), conf.kernel=c("linear","gaussian"), bandwidth=NULL, c1=NULL, c2=NULL, n.folds=5, maxiter=50, delta=10^-8, filter=T, minmaxsep=1e-3, n.cores=4, parallel=F)

cv.acSSPCA <- function(X,Y,A,lambda,npc,resp.kernel=c("linear","delta"),
                       conf.kernel=c("linear","gaussian"), bandwidth=NULL,
                       c1=NULL, c2=NULL, n.folds=5, maxiter=50, delta=10^-8,
                       filter=T, minmaxsep=1e-3, n.cores=4, parallel=F){
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
  }
  #if A is a vector, change it to a matrix
  if (is.null(dim(A))){
    A <- matrix(A, ncol=1)
  }
  
  n <- dim(X)[[1]]
  p <- dim(X)[[2]]
  
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
    c2.low <- 1
    c2.high <- sqrt(p) - minmaxsep
    c2.diff <- (c2.high-cd.low)/24
    c2 <- seq(from=c2.low,to=c2.high,by=c2.diff)
  }
  
  if(length(c2) < 2) stop("Must provide at least 2 sparse values")
  if(min(c2) < 1) stop("Cannot have minimum value less than 1")
  if(max(c2) > sqrt(p)) stop("Cannot have maximum value greater than the square root of the number of predictors")
  n.c2 <- length(c2)
  
  if(is.null(colnames(X))) colnames(X) <- paste0("X",1:p)
  if(is.null(colnames(Y))) colnames(Y) <- paste0("Y",1:dim(Y)[[2]])
  if(is.null(colnames(A))) colnames(A) <- paste0("A",1:dim(A)[[2]])
  
  x.names <- colnames(X)
  y.names <- colnames(Y)
  a.names <- colnames(A)
  
  df <- data.frame(Y,X,A)
  df.partition <- groupdata2::fold(data=df,k=n.folds)
  
  fold.arg <- c(1:n.folds)
  param.grid <- expand.grid(fold.arg,c2)
  colnames(param.grid) <- c("fold.arg","sp.arg")
  
  if(parallel){ 
    if(is.null(n.cores)) n.cores <- parallel::detectCores() 
    clust <- parallel::makeCluster(n.cores)
    metric.vec <- parallel::parApply(cl=clust,X=as.matrix(param.grid),1,
                                     cv.partition.acSSPCA,
                                     df.partition=df.partition,npc=npc,
                                     n.folds=n.folds,sparsity.type=sparsity.type,
                                     nonzero.loadings=NULL,
                                     sumabsv=NULL,kernel=kernel,
                                     niter=niter,trace=trace)
    
  } else{
    metric.vec <- apply(X=as.matrix(param.grid),1,cv.partition.acSSPCA,
                        df.partition=df.partition,npc=npc,
                        n.folds=n.folds,sparsity.type=sparsity.type,
                        nonzero.loadings=NULL,
                        sumabsv=NULL,kernel=kernel,niter=niter,trace=trace)
  }
  
  
  param.grid <- cbind(param.grid,metric.vec)
  metric.matrix <- mat.fill(param.grid=param.grid,n.sp=n.c2,n.folds=n.folds)
  
  cv.metric <- apply(metric.matrix,1,mean)
  cv.df <- data.frame(sparsity=c2,cv.metric=cv.metric)
  best.metric <- min(cv.df$cv.metric)
  best.sparse.param <- cv.df$sparsity[which(cv.df$cv.metric==best.metric)]
  
  return(list(cv.matrix=metric.matrix, cv.metrics=cv.df, best.metric=best.metric,best.sparse.param=best.sparse.param))
  
}
