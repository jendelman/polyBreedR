#' Extract read counts from AD string
#' 
#' Extract read counts from AD string
#' 
#' Only valid for a single ALT allele. 
#' 
#' @param AD array of AD strings
#' @param ALT TRUE or FALSE (= REF)
#' @param n.core number of cores
#' 
#' @return integer data with same dimensions as AD
#' @export
#' @importFrom parallel makeCluster stopCluster parApply clusterExport

ADsplit <- function(AD, ALT, n.core=1) {
  
  stopifnot(inherits(ALT,"logical"))
  
  f <- function(AD,k) {
    x <- strsplit(AD,split=",",fixed=T)
    m <- length(x)
    out <- integer(m)
    ok <- which(sapply(x,length)==2)
    out[ok] <- as.integer(sapply(x[ok],"[[",k))
    missing <- which(is.na(out))
    if (length(missing)>0)
      out[missing] <- 0
    return(out)
  }
  
  k <- as.integer(ALT)+1
  
  if (inherits(AD,"matrix")) {
    if (n.core > 1) {
      cl <- makeCluster(n.core)
      clusterExport(cl=cl,varlist=NULL)
      out <- parApply(cl=cl,X=AD,MARGIN=2,FUN=f,k=k)
      stopCluster(cl)
    } else {
      out <- apply(AD,2,f,k=k)
    }
    dimnames(out) <- dimnames(AD)
  } else {
    out <- f(AD,k)
    names(out) <- names(AD)
  }
  return(out)
}
