#' L2 norm calculation
#' @param vector vector to calculate l2 norm
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples getnorm2(vector)

getnorm2 <- function(vector){
  return(sqrt(sum(vector^2)))
}