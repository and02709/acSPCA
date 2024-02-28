#' Generates response kernel.  Can be either linear or delta kernel
#' @param Y Response vector
#' @param n number of observations in response vector Y
#' @param kernel Type of kernel.  Either linear or delta
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples respkernel(Y, n, kernel)

respkernel <- function(Y, n, kernel){
  
  if (kernel=="linear"){
    L <- Tcrossprod(Y,Y)
  } else if (kernel=="delta"){
    yf <- as.factor(Y)
    L<-matrix(0,n,n)
    for (i in levels(yf)){
      tmp<-yf==i
      L[tmp,tmp]<-1
    } 
  } else {
    stop("Please select a valid kernel, linear kernel or delta kernel")
  } 
  return(L)
}