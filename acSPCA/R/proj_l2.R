#' Internal function to perform l2n penalization on vector
#' @param v vector to under go l2 norm normalization
#' @param c constant to multiply l2n norm 1 vector by
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples proj_l2function(v, c)

proj_l2 <- function(v, c){
  return(v/getnorm2(v)*c)    
}