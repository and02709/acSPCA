#' Internal function to perform L1 penalization on vector
#' @param M Confounding matrix product with predictors 
#' @param S Singular values of M matrix 
#' @param V right vectors of M matrix 
#' @param utx crossproduct of u vector and Psi matrix
#' @param c2 sum of absolute value loadings for desired sparse vector 
#' @param v_ini initial vector to undergo L1 penalization and projection 
#' @param tmin iterative lower bound 
#' @param tmax iterative upper bound 
#' @param iternum number of iterations
#' @param minmaxsep determines the acceptable algorithm termination acceptable different
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples proj( M, S, V, utx, c1, c2, v_ini, tmin, tmax, iternum, minmaxsep)


proj <- function( M, S, V, utx, c1, c2, v_ini, tmin, tmax, iternum, minmaxsep){
  epsilon <- minmaxsep/iternum; epsilon1 <- 10^-3/iternum; epsilon2 <- 10^-4/iternum; maxiter <- round(50*sqrt(iternum))
  tlow <- tmin; tup <- tmax; v <- v_ini
  while (tup - tlow >= tup*epsilon){
    t <- (tup + tlow)/2
    iter <- 1; converge <- 0; 
    while (iter<=maxiter && converge==0){
      #cat(iternum, tlow, tup, t, iter,"\n")
      #project to l1 ball
      if (sum(abs(v)) > c2){
        v <- proj_l1(v, c2);  
      } 
      #project to l2 ball
      if (getnorm2(v) > 1){
        v <- proj_l2( v, 1);  
      }
      #project to hyperplane
      if (t - utx%*%v > 0){
        v <- proj_hp( v, t(utx), t)
      }
      #project to elliptic
      P <- crossprod(V, v)^2; 
      if (getnorm2(M%*%v) > sqrt(c1)){
        l <- proj_l2g_getl( S, P, c1)
        ####when there is only one non-zero singular value
        if (length(S)==1){
          v <- v - V*(1/(1 + 1/(l*S^2)))*as.numeric(crossprod(V,v))
        } else {
          v <- v - V%*%(diag(1/(1 + 1/(l*S^2)))%*%crossprod(V,v))  
        }
        
      }
      #check constrain
      if (t - utx%*%v <= t*epsilon & sum(abs(v)) <= c2*(1 + epsilon1) & getnorm2(v) <= 1 + epsilon2){
        converge <- 1  
      }
      iter <- iter + 1
    }
    
    if (converge == 1){
      tlow <- utx%*%v
      vcon <- v
    } else {
      tup <- t  
    }
  }
  return(list(v = vcon, t = utx%*%vcon, conv = converge))
}