#' Calculates the kernel for the confounding variables
#' @param Y Confounding variable matrix
#' @param kernel Type of kernel.  Either linear or gaussian 
#' @param bandwidth scaling parameter for gaussian kernel
#' @param scaleY scales the confounding variables
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples respkernel(Y, n, kernel)

calkernel <- function(Y, kernel, bandwidth, scaleY=F){
  Y <- scale(Y, center = F, scale = scaleY)  
  ####missing data
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  
  if (kernel=="linear"){
    K <- tcrossprod(Y)
  } else if (kernel=="gaussian"){
    if (is.null(bandwidth)==T){
      stop("For gaussian kernel, please specify the bandwidth") 
    } else{
      K <- as.matrix(dist(Y, method = "euclidean"))
      K <- exp(-K^2/2/bandwidth^2)  
    }
  } else {
    stop("Please select a valid kernel, linear kernel or gaussian kernel")
  } 
  return(K)
}