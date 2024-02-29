#' Internal function to perform l1 penalization
#' @param v vector to under go l1 penalty
#' @param c l1 penalty term
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples proj_l1(v, c)

proj_l1 <- function(v, c){
  if (sum(abs(v)) < c){
    return(v)
  }
  u <- sort(abs(v), decreasing=T)
  sv <- cumsum(u)
  rho <- max(which(u > (sv - c) / (1:length(u))))
  theta <- max(0, (sv[rho] - c) / rho)
  return( sign(v) * pmax(abs(v) - theta, 0) )
}