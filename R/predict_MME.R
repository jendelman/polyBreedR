#' Compute BLUPs by solving the Mixed Model Equations
#' 
#' Compute BLUPs by solving the Mixed Model Equations
#' 
#' BLUPs are computed at the average value of the fixed effects. If \code{weights} is used, the names must exactly match the names of the kernels in \code{data}. Using the argument \code{exclude.id}, a subset of the population can be excluded, to enable cross-validation. The function \code{\link{Stage2}} can be used to create a suitable object of class \code{\link{MME}}.
#' 
#' @param data Variable of class \code{\link{MME}}
#' @param weights Named vector of weights for the genetic effects in BLUP. Default is 1 for all effects.
#' @param exclude.id Vector of individuals to exclude from the training set (optional)
#' 
#' @return data frame with columns id,blup,r2
#' @import Matrix
#' @export
#' 
predict_MME <- function(data,weights=NULL,exclude.id=NULL) {
  
  stopifnot(!is.na(data@y))
  n.random <- length(data@kernels)
  stopifnot(n.random > 0)
  if (is.null(weights)) {
    weights <- rep(1,n.random)
    names(weights) <- names(data@kernels)
  }
  stopifnot(any(weights > 0))
  stopifnot(n.random==length(weights))
  stopifnot(setequal(names(weights),names(data@kernels)))
  
  weights <- Matrix(weights[match(names(weights),names(data@kernels))],nrow=1,ncol=n.random)
  idK <- lapply(data@kernels,rownames)
  id <- colnames(data@Z)
  n <- length(id)
  stopifnot(sapply(idK,setequal,y=id))
  
  if (!is.null(exclude.id)) {
    exclude.id <- unique(as.character(exclude.id))
    ix <- which(is.element(id,exclude.id))
    if (length(ix) > 0) {
      remove <- as.integer(unlist(apply(as.matrix(data@Z[,ix]),2,function(x){which(x==1)})))
      if (length(remove)>0) {
        data@y <- data@y[-remove]
        data@X <- data@X[-remove,]
        j <- which(apply(data@X,2,sum)==0)
        if (length(j)>0) {
          warning("Excluded set eliminated an environment")
          data@X <- data@X[,-j]
        }
        data@Z <- data@Z[-remove,]
        data@Rmat <- data@Rmat[-remove,-remove]
      }
    }
  }
  
  Z <- kronecker(Matrix(1,nrow=1,ncol=n.random),data@Z)
  RX <- solve(data@Rmat,data@X)
  RZ <- solve(data@Rmat,Z)
  Ry <- solve(data@Rmat,data@y)
  Q11 <- crossprod(data@X,RX)
  Q12 <- crossprod(data@X,RZ)
  Q22 <- crossprod(Z,RZ)+direct_sum(lapply(data@kernels,solve))
  Q <- rbind(cbind(Q11,Q12),cbind(t(Q12),Q22))
  b1 <- as.numeric(crossprod(data@X,Ry))
  b2 <- as.numeric(crossprod(Z,Ry))
  b <- c(b1,b2)
  C <- solve(Q)
  ans <- as.numeric(C%*%b)
  n.fix <- ncol(data@X)
  fixed <- sum(apply(data@X,2,mean)*ans[1:n.fix])
  M <- kronecker(weights,Diagonal(n))
  numer <- diag(M%*%C[-(1:n.fix),-(1:n.fix)]%*%t(M))
  denom <- diag(M%*%direct_sum(data@kernels)%*%t(M))
  return(data.frame(id=id,blup=as.numeric(M%*%ans[-(1:n.fix)])+fixed,r2=1-numer/denom,stringsAsFactors = F))
}
