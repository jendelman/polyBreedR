#' Predict allele dosage
#' 
#' Returns allele dosage as the maximum a posteriori (MAP) genotype based on a normal mixture model (NMM)
#' 
#' The input values should fall within the interval [0,1]. The vector of model parameters should be ordered as follows: means, standard deviations, mixture probabilities. If the posterior probability of the MAP genotype is less than `min.posterior`, then `NA` is returned for that sample. By default, an arcsin square root transformation is applied to the input values, as this is the approach used by the `fitPoly` package. 
#' 
#' @param x vector of input values 
#' @param params vector of parameters for the normal mixture model
#' @param ploidy 2 or 4
#' @param min.posterior minimum posterior probability for a genotype call to be made
#' @param transform TRUE/FALSE whether to apply an arcsin square root transformation
#' 
#' @keywords internal
#' @return vector of allele dosages (0,1,2,..ploidy)
#' 
#' @importFrom stats dnorm

predict_NMM <- function(x,params,ploidy,min.posterior=0,transform=TRUE) {
  f <- function(x) {
    z <- ifelse(x < 0,0,x)
    return(ifelse(z>1,1,z))
  }
  n <- length(x)
  x <- f(x)
  if (transform) {x <- asin(sqrt(x))}
  
  likelihood <- matrix(0,nrow=n,ncol=ploidy+1)
  for (k in 0:ploidy) {
    likelihood[,k+1] <- dnorm(x,mean=params[k+1],sd=params[(ploidy+1)+(k+1)])
  }
  mixing.probs <- params[(ploidy+1)*2+(0:ploidy)+1]
  mixing.probs <- mixing.probs/sum(mixing.probs) #normalize
  post <- likelihood * kronecker(matrix(1,nrow=n,ncol=1),matrix(mixing.probs,nrow=1,ncol=(ploidy+1)))
  MAP.geno <- apply(post,1,function(p){
    p <- p/sum(p)
    maxp <- max(p)
    if (is.na(maxp)) {return(NA)}
    if (maxp>=min.posterior) {
      return(which.max(p)-1)
    } else {
      return(NA)
    }
  })
  return(MAP.geno)
}
