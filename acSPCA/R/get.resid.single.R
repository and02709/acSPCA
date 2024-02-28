#' For loop residualization
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param a nxq confounding matrix
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples get.resid.single(X,Y,a)

get.resid.single = function(X, Y, a){
  X.tilde = c()
  X.star = c()
  for (j in 1:ncol(X)){
    lm.X.j = lm(X[,j] ~ Y + a)
    X.tilde = cbind(X.tilde, (X[,j] - Y%*%lm.X.j$coefficients[startsWith(names(lm.X.j$coefficients),"Y")]))
    X.star = cbind(X.star, (X[,j] - a%*%lm.X.j$coefficients[startsWith(names(lm.X.j$coefficients),"a")]))
  }
  
  return(list(X.tilde = X.tilde,
              X.star = X.star))
}