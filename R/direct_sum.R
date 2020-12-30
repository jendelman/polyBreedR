#' Direct Sum
#' 
#' Direct Sum
#' 
#' Computes the direct sum of the matrices in \code{x}
#' 
#' @param x list of matrices
#' 
#' @return Sparse Matrix
#' @import Matrix
#' 
#' @keywords internal
#' 
direct_sum <- function(x) {
  n <- length(x) 
  m <- sapply(x,nrow)
  m.cumulative <- apply(array(1:n),1,function(k){sum(m[1:k])})
  i=1
  z <- expand.grid(col=1:m[i],row=1:m[i])
  out <- data.frame(row=z$row,col=z$col,value=as.vector(x[[i]]))
  for (i in 2:n) {
    z <- expand.grid(col=(1:m[i])+m.cumulative[i-1],row=(1:m[i])+m.cumulative[i-1])
    tmp <- data.frame(row=z$row,col=z$col,value=as.vector(x[[i]]))
    out <- rbind(out,tmp)
  }
  #return(out)
  return(sparseMatrix(i=out[,1],j=out[,2],x=out[,3]))
}
