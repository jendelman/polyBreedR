#' Impute missing marker data 
#' 
#' Impute marker data based on the population mean or mode
#' 
#' Missing values are imputed with either the population mean or mode (most frequent value) for each marker
#' 
#' @param geno Matrix of allele dosages with dimensions markers x indiv 
#' @param method Either "mean" or "mode"
#' 
#' @return Imputed genotype matrix (markers x indiv)
#' 
#' @export

impute <- function(geno,method) {
  stopifnot(method %in% c("mode","mean"))
  
  impute.mode <- function(x) {
    miss <- which(is.na(x))
    if (length(miss)>0) {
      x[miss] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  impute.mean <- function(x) {
    miss <- which(is.na(x))
    if (length(miss)>0) {
      x[miss] <- mean(x[-miss])
    }
    return(x)
  }
  if (method=="mode") {
    return(t(apply(geno,1,impute.mode)))
  }
  if (method=="mean") {
    return(t(apply(geno,1,impute.mean)))
  }
}
