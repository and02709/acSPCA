#' CV function for determining sparsity penalty c2
#' @param arg.sparse vector for folds and sparse parameter
#' @param df.partition list of folds
#' @param lambda penalty for adjusting for confounding 
#' @param npc number of principal components 
#' @param n.folds number of folds for cross-validation
#' @param c1 norm of quadratic predictor and confounding matrix product
#' @param resp.kernel kernel of response vector which can be linear or delta 
#' @param conf.kernel kernel of confounding matrix which can be linear or gaussian 
#' @param bandwidth gaussian kernel argument 
#' @param maxiter maximum number of iterations for the algorithm 
#' @param delta offset parameter
#' @param filter parameter to determine filtering
#' @param minmaxsep determines the acceptable algorithm termination acceptable different 
#' @param x.names column names for predictors 
#' @param y.names name of response variable 
#' @param a.names name of confounding variables
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples cv.acSSPCA <- function(X,Y,A,lambda,npc,resp.kernel=c("linear","delta"), conf.kernel=c("linear","gaussian"), bandwidth=NULL, c1=NULL, c2=NULL, n.folds=5, maxiter=50, delta=10^-8, filter=T, minmaxsep=1e-3, n.cores=4, parallel=F)

cv.partition.acSSPCA <- function(arg.sparse, df.partition, lambda, npc, n.folds, 
                                 resp.kernel=c("linear","delta"),
                                 conf.kernel=c("linear","gaussian"),
                                 bandwidth=NULL, c1=NULL, maxiter=50, 
                                 delta=10^-8, filter=T, minmaxsep=1e-3,
                                 x.names=NULL, y.names=NULL, a.names=NULL){
  test.index <- arg.sparse[[1]]
  c2 <- arg.sparse[[2]]
  
  X <- df.partition |> dplyr::filter(.folds!=test.index) |> dplyr::ungroup() |> dplyr::select(tidyselect::all_of(x.names)) |> as.matrix()
  xtest <- df.partition |> dplyr::filter(.folds==test.index) |> dplyr::ungroup() |> dplyr::select(tidyselect::all_of(x.names)) |> as.matrix()
  Y <- df.partition |> dplyr::filter(.folds!=test.index) |> dplyr::ungroup() |> dplyr::select(tidyselect::all_of(y.names)) |> as.matrix()
  ytest <- df.partition |> dplyr::filter(.folds==test.index) |> dplyr::ungroup() |> dplyr::select(tidyselect::all_of(y.names)) |> as.matrix()
  A <- df.partition |> dplyr::filter(.folds!=test.index) |> dplyr::ungroup() |> dplyr::select(tidyselect::all_of(a.names)) |> as.matrix()
  atest <- df.partition |> dplyr::filter(.folds==test.index) |> dplyr::ungroup() |> dplyr::select(tidyselect::all_of(a.names)) |> as.matrix()
  
  xmeans <- colMeans(X)
  xtrain <- t(apply(X, 1, function(x) x-xmeans))
  xtest <- t(apply(xtest, 1, function(x) x-xmeans))
  
  acSSPCA.obj <- acSPCA::acSSPCAmulti(X=X,Y=Y,A=A,lambda=lambda,npc=npc,resp.kernel=resp.kernel, conf.kernel=conf.kernel, bandwidth=bandwidth, 
                                      c1=c1, c2=c2, maxiter=maxiter, delta=delta, filter=filter, minmaxsep=minmaxsep)
  
  V <- acSSPCA.obj$vectors
  
  mod <- acSPCA::model.build.SSPCA(Y=Y,X=X,V=V,kernel=resp.kernel)
  metric <- acSPCA::model.metric.SSPCA(Y=ytest,X=xtest,V=V,mod=mod,kernel=resp.kernel)
  return(metric)
}