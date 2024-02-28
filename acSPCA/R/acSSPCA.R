#' Sparse Supervised Adjusted for Confounding Principal Component Analysis
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param A nxq confounding variable matrix
#' @param X4A alternative confounding matrix starting point
#' @param x4Y alternative response variable starting point 
#' @param c1 norm of quadratic predictor and confounding matrix product
#' @param c2 sum of absolute value loadings for desired sparse vector 
#' @param v_ini non-sparse initial vector 
#' @param resp.kernel kernel of response vector which can be linear or delta 
#' @param conf.kernel kernel of confounding matrix which can be linear or gaussian 
#' @param bandwidth gaussian kernel argument 
#' @param maxiter maximum number of iterations for the algorithm 
#' @param delta offset parameter
#' @param filter parameter to determine filtering
#' @param minmaxsep determines the acceptable algorithm termination acceptable different
#' @keywords Sparse Adjusted Confounding Supervised Principal Component Analysis
#' @export
#' @examples acSSPCA(X, Y, A, X4A=NULL, X4Y=NULL, c1=NULL, c2, v_ini, resp.kernel=c("linear","delta"), conf.kernel=c("linear", "gaussian"), bandwidth=NULL, maxiter=50, delta=10^-8, filter=T, minmaxsep=1e-3)

acSSPCA <- function(X, Y, A, X4A=NULL, X4Y=NULL, c1=NULL, c2, 
                    v_ini, resp.kernel=c("linear","delta"), 
                    conf.kernel=c("linear", "gaussian"), bandwidth=NULL, 
                    maxiter=50, delta=10^-8, filter=T, minmaxsep=1e-3){
  
#if Y is a vector, change it to a matrix
if (is.null(dim(Y))){
  Y <- matrix(Y, ncol=1)
}
#if A is a vector, change it to a matrix
if (is.null(dim(A))){
  A <- matrix(A, ncol=1)
}
  
# Determine number of observations
n <- dim(Y)[[1]]

# Determine number of predictors
p <- dim(X)[[2]]

#check whether a whole row in X is missing
Xmis <- apply(X, 1, function(row){sum(!is.na(row))})
if (sum(Xmis==0)){
  stop(paste("The following samples in X is missing, please remove them in X and Y: rows ", paste(which(Xmis==0), collapse =  " "), sep=""))
}

#check whether the number of samples in X, Y, and A match
if (dim(X)[1]!=n){
  stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
}

if (dim(X)[1]!=dim(A)[1]){
  stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
}

if (n!=dim(A)[1]){
  stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
}

####check whether the number of samples in X4Y and Y match
if (!is.null(X4Y)){
  if (dim(X4Y)[1]!=n){
    stop("The numbers of samples in X4Y ( nrow(X4Y) ) and Y ( nrow(Y) ) do not match")
  }  
}
#scale data
# X <- scale(X, center = centerX, scale = scaleX)  
# Y <- scale(Y, center = centerY, scale = scaleY) 
# A <- scale(A, center = centerA, scale = scaleA)

#missing data
X[is.na(X)] <- mean(X, na.rm=T)
Y[is.na(Y)] <- mean(Y, na.rm=T)
A[is.na(A)] <- mean(A, na.rm=T)
if (!is.null(X4A)){
  X4A[is.na(X4A)] <- mean(X4A, na.rm=T)
}
if (!is.null(X4Y)){
  X4Y[is.na(X4Y)] <- mean(X4Y, na.rm=T)
}

Ka <- calkernel(A, conf.kernel, bandwidth)
Ky <- respkernel(Y, n, resp.kernel)
H <- diag(x=1, nrow=n) - 1/n*rep(1,n)%*%t(rep(1,n))

if (is.null(c1)){
  if (is.null(X4A)){
    HXv <- H%*%X%*%as.matrix(v_ini)
    c1 <- as.numeric(crossprod(HXv, Ka%*%HXv))  
  } else {
    HXv <- H%*%X4Y%*%as.matrix(v_ini)
    c1 <- as.numeric(crossprod(HXv, Ka%*%HXv)) 
  }
}
#check if either c1 or c2 is 0
if (c1==0 | c2==0){
  return(list(v=as.matrix(v_ini*0), u=X%*%as.matrix(v_ini*0), converge=1))  
}

#Decomposition on the kernel matrix to get the M matrix and Psi matrix
M <- get.psi(X=X,L=Ka,n=n,p=p)
Psi <- get.psi(X=X,L=Ky,n=n,p=p)

####svd on M
tmp = svd(M);
S <- tmp$d; V <- tmp$v

####initialization
v <- v_ini;
valA <- c();
devia <- c();
iter <- 1; converge <- 0;
val <- -Inf;
####iteration
while (converge==0 & iter<= maxiter){
  #cat(iter," ")
  ####update u
  tmp <- Psi%*%v
  u <- tmp/getnorm2(tmp)
  ####update v
  utx <- crossprod(u,Psi)
  tmin <- 0; tmax <- getnorm2(utx)
  tmpproj <- proj( M, S, V, utx, c1, c2, v, tmin, tmax, iter, minmaxsep)
  vnew <- tmpproj$v
  #check convergence
  devia <- c(devia, getnorm2(v-vnew))
  valnew <- tmpproj$t
  if (valnew - val <= delta){
    converge <- 1
  }
  iter <- iter + 1
  v <- vnew
  val <- valnew
  valA <- c(valA, val)    
}
####Because of the bioconvexity of the problem and numerical issues, sometimes there are a lot of small entries in v very close to 0
if (filter==T){
  v <- as.matrix(topthresfilter(v, top=0.01, alpha=0.5*10^-4))
}
return(list(v=v, u=u, converge=converge))
}