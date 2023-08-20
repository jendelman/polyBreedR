#' Convert GT to DS
#' 
#' Convert GT to DS
#' 
#' @param GT GT string
#' 
#' @return numeric DS
#' @export

GT2DS <- function(GT) {
  #check for missing data
  options(warn = -1)
  geno <- sapply(strsplit(as.character(GT),split="/",fixed=T),function(x){sum(as.integer(x))})
  options(warn = 0)
  return(geno)
}
