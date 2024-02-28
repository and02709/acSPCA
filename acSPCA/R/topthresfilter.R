#' This function converts very low magnitude loadings into zero
#' @param v vector to be converted
#' @param top theshold value
#' @param alpha theshold limit
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples topthresfilter(v, top, alpha)

topthresfilter <- function(v, top, alpha){
  thres <- mean(sort(abs(v), decreasing = T)[1:round(length(v)*top)])*alpha
  v[abs(v)<=thres] <- 0
  return(v)
}