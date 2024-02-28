#' Extract single vector for SIMPLS algorithm
#' @param X nxp predictor matrix
#' @param H centering matrix 
#' @param Ky response kernel matrix
#' @param Ka confounding kernel matrix
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples single.vector.decomp <- function(X,H,Ky,Ka,lambda)

single.vector.decomp <- function(X,H,Ky,Ka,lambda){
  return(eigs_sym(mat.mult(mat.mult(mat.mult(Crossprod(X,H),Ky),H),X) - lambda*mat.mult(mat.mult(mat.mult(Crossprod(X,H),Ka),H),X), k=1, which="LA"))
}