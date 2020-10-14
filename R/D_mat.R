#' Dominance genomic relationships
#' 
#' Coefficients and relationship matrix for digenic dominance effects with bi-allelic markers
#' 
#' Digenic dominance effects are based on the traditional orthogonal decomposition of genetic variance in panmictic populations (Fisher 1918; Kempthorne 1957; Endelman et al. 2018). The D matrix is computed from the coefficients and scaling factor according to D = tcrossprod(coeff/scale). Missing genotype data is replaced with the population mean.
#' 
#' @references Fisher (1918) Trans. Roy. Soc. Edin. 52:399-433.
#' @references Kempthorne (1957) An Introduction to Genetic Statistics.
#' @references Endelman et al. (2018) Genetics 209:77-87.
#'
#' @param geno Matrix of allele dosages: markers x indiv
#' @param ploidy 2 or 4
#' 
#' @return List containing
#' \describe{
#' \item{coeff}{Coefficients of the marker effects (dim: indiv x marker)}
#' \item{scale}{Scaling factor between markers and indiv}
#' \item{mat}{D matrix}
#' }
#' 
#' @export

D_mat <- function(geno,ploidy) {
  stopifnot(ploidy %in% c(2,4))
  m <- nrow(geno)
  n <- ncol(geno)
  geno <- impute(geno,method="mean")
  
  p <- matrix(apply(geno,1,mean)/ploidy,nrow=m,ncol=n,byrow=F)
  if (ploidy==4) {
    coeff <- -12*p^2+(6*p+1)*geno-geno^2
    scale <- sqrt(sum(24*p[,1]^2*(1-p[,1])^2))
  } else {
    #ploidy=2
    coeff <- -2*p^2+(2*p+1)*geno-geno^2
    scale <- sqrt(sum(4*p[,1]^2*(1-p[,1])^2))
  }
  return(list(coeff=t(coeff),scale=scale,mat=crossprod(coeff/scale)))
}