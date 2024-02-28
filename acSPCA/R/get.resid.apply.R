#' Uses apply function to perform residual calculation
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param a nxq confounding matrix
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples get.resid.apply(X,Y,a)

get.resid.apply <- function(X, Y, a){
  temp <- apply(X=X,2,res.calc,Y=Y,a=a)
  X.star <- unlist(lapply(X=temp, function(x) return(x$star)))
  X.tilde <- unlist(lapply(X=temp, function(x) return(x$tilde)))
  return(list(X.tilde=X.tilde,X.star=X.star))
}