#' Additive genomic relationships
#' 
#' Relationship matrix for additive effects with bi-allelic markers
#' 
#' Additive effects are based on the traditional orthogonal decomposition of genetic variance in panmictic populations (Fisher 1918; Kempthorne 1957; Endelman et al. 2018). Missing genotype data is replaced with the population mean. Equivalent to VanRaden Method 1 for diploids.
#' 
#' If \code{p.ref} is NULL, the current population is used as the reference population.
#' 
#' @references Fisher (1918) Trans. Roy. Soc. Edin. 52:399-433.
#' @references Kempthorne (1957) An Introduction to Genetic Statistics.
#' @references VanRaden (2008) J. Dairy Sci 91:4414-4423.
#' @references Endelman et al. (2018) Genetics 209:77-87.
#'  
#' @param geno Matrix of allele dosages (markers x indiv)
#' @param ploidy Any even integer (2,4,6,...)
#' @param p.ref reference population frequency
#' 
#' @return G matrix 
#' 
#' @export

G_mat <- function(geno, ploidy, p.ref=NULL) {
  m <- nrow(geno)
  n <- ncol(geno)
  if (is.null(p.ref)) {
    p.ref <- apply(geno,1,mean,na.rm=T)/ploidy
  }
  Pmat <- kronecker(matrix(p.ref,nrow=1,ncol=m),matrix(1,nrow=n,ncol=1))
  coeff <- t(geno) - ploidy * Pmat 
  coeff[is.na(coeff)] <- 0
  dimnames(coeff) <- dimnames(geno)[c(2,1)]
  scale <- ploidy*sum(p.ref*(1-p.ref))
  return(tcrossprod(coeff)/scale)
}
