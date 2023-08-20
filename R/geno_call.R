#' Genotype calls 
#' 
#' Genotype calls based on a normal mixture model
#' 
#' The first column of the CSV input file should be the SNP ID, followed by columns for the normal distribution means, standard deviations, and mixture probabilities. Genotype calls are based on the maximum a posteriori (MAP) method. If the posterior probability of the MAP genotype is less than \code{min.posterior}, then \code{NA} is returned for that sample. By default, an arcsin square root transformation is applied to the input values to match the approach used by R package fitPoly. To use a tetraploid mixture model for diploid samples, set \code{sample.ploidy = 2} and \code{model.ploidy = 4}. 
#' 
#' @param data matrix (markers x id) of input values for the normal mixture model
#' @param filename CSV filename with the model parameters
#' @param model.ploidy 2 or 4 (default)
#' @param sample.ploidy 2 or 4 (default)
#' @param min.posterior minimum posterior probability (default 0) for genotype call
#' @param transform TRUE (default) or FALSE whether to apply arcsin square root transformation
#' 
#' @return matrix of allele dosages (0,1,2,..ploidy) with dimensions markers x individuals
#' 
#' @export
#' @importFrom utils read.csv

geno_call <- function(data,filename,model.ploidy=4L,sample.ploidy=4L,min.posterior=0,transform=TRUE) {
  
  model.params <- as.matrix(read.csv(file=filename,as.is=T,row.names = 1))
  if (is.null(model.ploidy)) {
    model.ploidy <- ncol(model.params)/3 - 1
  }
  if (is.null(sample.ploidy)) {
    sample.ploidy <- model.ploidy
  }
  
  stopifnot(sample.ploidy <= model.ploidy)
  stopifnot(sample.ploidy %in% c(2L,4L))
  stopifnot(model.ploidy %in% c(2L,4L))

  iv <- which(rownames(data) %in% rownames(model.params))
  snp.id <- rownames(data)[iv]
  sample.id <- colnames(data)
  m <- length(iv)
  data <- lapply(1:m,function(i){data[iv[i],]}) #make it a list
  ix <- match(snp.id,rownames(model.params))
  if (sample.ploidy == model.ploidy) {
    model.params <- lapply(1:m,function(i){model.params[ix[i],]})
  } else {
    #2x sample with 4x model
    model.params <- lapply(1:m,function(i){model.params[ix[i],c(1,3,5,6,8,10,11,13,15)]})
  }
  
  geno <- mapply(FUN=predict_NMM,x=data,params=model.params,ploidy=sample.ploidy,min.posterior=min.posterior,transform=transform)
  if (length(sample.id)==1) {
    geno <- matrix(geno,ncol=1)
  } else {
    geno <- t(geno)
  }
  colnames(geno) <- sample.id
  rownames(geno) <- snp.id
  return(geno)
}
