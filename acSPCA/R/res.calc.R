#' Residual calculation
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param a nxq confounding matrix
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples res.calc(X,Y,a)

res.calc <- function(x, Y, a){
  lm.X.j <- lm(x ~ Y + a)
  X.tilde <- (x - Y%*%lm.X.j$coefficients[startsWith(names(lm.X.j$coefficients),"Y")])
  X.star <- (x - a%*%lm.X.j$coefficients[startsWith(names(lm.X.j$coefficients),"a")])
  return(list(star=X.star, tilde=X.tilde))
}