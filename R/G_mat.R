#' Additive genomic relationships
#' 
#' Relationship matrix for additive effects with bi-allelic markers
#' 
#' Additive effects are based on the traditional orthogonal decomposition of genetic variance in panmictic populations (Fisher 1918; Kempthorne 1957; Endelman et al. 2018). Missing genotype data is replaced with the population mean. 
#' 
#' @references Fisher (1918) Trans. Roy. Soc. Edin. 52:399-433.
#' @references Kempthorne (1957) An Introduction to Genetic Statistics.
#' @references Endelman et al. (2018) Genetics 209:77-87.
#'  
#' @param geno Matrix of allele dosages (markers x indiv)
#' @param ploidy Any even integer (2,4,6,...)
#' 
#' @return G matrix 
#' 
#' @export

G_mat <- function(geno,ploidy) {
  m <- nrow(geno)
  n <- ncol(geno)
  coeff <- scale(t(geno),scale=F)[1:n,1:m]
  coeff[is.na(coeff)] <- 0
  p <- apply(geno,1,mean,na.rm=T)/ploidy
  scale <- sqrt(sum(ploidy*p*(1-p)))
  return(tcrossprod(coeff/scale))
}
