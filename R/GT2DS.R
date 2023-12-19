#' Convert GT to ALT allele dosage (DS)
#' 
#' Convert GT to ALT allele dosage (DS)
#' 
#' If \code{diploidize} is TRUE, data are recoded as a diploid {0,1,2}.
#' 
#' @param GT GT string
#' @param diploidize TRUE/FALSE
#' @param n.core number of cores
#' 
#' @return integer data with same dimensions as GT
#' @export
#' @importFrom parallel makeCluster stopCluster parApply clusterExport

GT2DS <- function(GT, diploidize=FALSE, n.core=1) {
  
  G2D <- function(GT, diploidize) {
    missing <- which(GT==".")
    ploidy <- sapply(gregexpr("/",GT,fixed=T),length)+1
    if (length(missing) > 0)
      ploidy[missing] <- NA
    options(warn = -1)
    geno <- sapply(strsplit(as.character(GT),split="/",fixed=T),function(x){sum(as.integer(x))})
    options(warn = 0)
    if (diploidize) {
      geno <- ifelse(geno==0,0,ifelse(geno==ploidy,2,1))
    }
    return(geno)
  }
  
  if (inherits(GT,"matrix")) {
    if (n.core > 1) {
      cl <- makeCluster(n.core)
      clusterExport(cl=cl,varlist=NULL)
      geno <- parApply(cl=cl,X=GT,MARGIN=2,FUN=G2D,diploidize=diploidize)
      stopCluster(cl)
    } else {
      geno <- apply(GT,2,G2D,diploidize=diploidize)
    }
    dimnames(geno) <- dimnames(GT)
  } else {
    geno <- G2D(GT,diploidize=diploidize)
    names(geno) <- names(GT)
  }
  return(geno)
}
