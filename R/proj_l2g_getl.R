#' Projection onto hyperplane
#' @param S Singular values
#' @param P Projection matrix for v
#' @param c norm of crossprod between predictor matrix and confounding matrix
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples proj_l2g_getl( S, P, c)

proj_l2g_getl <- function( S, P, c){
  S2 <- S^2;
  epsilon <- 10^-10; maxiter <- 200;
  l <- 10; lnew <- 0
  iter <- 1;
  while (abs(l - lnew)>=epsilon && iter <= maxiter){
    l <- lnew  
    d1 <- sum(P*S2/((1 + l*S2)^2)) - c
    d2 <- -2*sum( P*S^4/ ( (1 + l*S2 )^3 ) )
    lnew <- l - d1/d2
  }
  return(lnew)  
}