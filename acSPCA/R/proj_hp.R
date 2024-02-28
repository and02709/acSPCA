#' Projection onto hyperplane
#' @param v vector to under go projection
#' @param w weight vector
#' @param t iterative value
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples proj_hp( v, w, t)

proj_hp <- function( v, w, t){
  t <- as.numeric(t)
  return( v - w*((sum(w*v) - t)/getnorm2(w)^2) )
}