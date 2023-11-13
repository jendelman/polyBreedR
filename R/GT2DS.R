#' Convert GT to DS
#' 
#' Convert GT to DS
#' 
#' @param GT GT string
#' @param diploidize TRUE/FALSE
#' 
#' @return numeric DS
#' @export

GT2DS <- function(GT,diploidize=FALSE) {
  #check for missing data
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
