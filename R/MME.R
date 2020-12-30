#' Solve the Mixed Model Equations
#' 
#' Solve the Mixed Model Equations
#' 
#' Details
#' 
#' @param y response (vector)
#' @param X design matrix for the fixed effects
#' @param Z list of design matrices for the random effects
#' @param Vu list of covariance matrices for the random effects
#' @param Rmat residual variance-covariance matrix
#' 
#' @return List containing
#' \describe{
#' \item{blue}{BLUE for the average fixed effect (numeric)}
#' \item{blup}{list of data frames, one for each random effect, with two variables: BLUP, r2}
#' }
#' 
#' @keywords internal
#' 
MME <- function(y,X,Z,Vu,Rmat) {
  stopifnot(!is.na(y))
  stopifnot(length(Z)==length(Vu))
  stopifnot(mapply(FUN=function(n1,n2){n1==n2},lapply(Z,ncol),lapply(Vu,ncol)))
  n.random <- length(Z)
  stopifnot(n.random > 0)
  Z2 <- NULL
  m <- integer(n.random+1)
  for (i in 1:n.random) {
    m[i+1] <- ncol(Z[[i]]) + m[i]
    Z2 <- cbind(Z2,Z[[i]])
  }
  RX <- solve(Rmat,X)
  RZ <- solve(Rmat,Z2)
  Ry <- solve(Rmat,y)
  Q11 <- crossprod(X,RX)
  Q12 <- crossprod(X,RZ)
  Q21 <- t(Q12)
  Q22 <- crossprod(Z2,RZ)+as.matrix(direct_sum(lapply(Vu,solve)))
  Q <- rbind(cbind(Q11,Q12),cbind(Q21,Q22))
  b1 <- crossprod(X,Ry)
  b2 <- crossprod(Z2,Ry)
  b <- c(b1,b2)
  C <- solve(Q)
  ans <- as.numeric(C%*%b)
  PEV <- diag(C)
  n.fix <- ncol(X)
  fixed <- sum(apply(X,2,mean)*ans[1:n.fix])
  BLUP <- vector("list",n.random)
  for (i in 1:n.random) {
    iv <- n.fix+(m[i]+1):m[i+1]
    BLUP[[i]] <- data.frame(BLUP=ans[iv],r2=1-PEV[iv]/diag(Vu[[i]]))
    rownames(BLUP[[i]]) <- colnames(Z[[i]])
  }
  return(list(blue=fixed,blup=BLUP))
}
