#' Genomic relationship matrix
#' 
#' Genomic relationship matrix 
#' 
#' For \code{method}="VR1", Method 1 of VanRaden (2008) is used, and its polyploid extension (Endelman et al. 2018). Missing data is replaced with the population mean for each marker. If \code{p.ref} is NULL, the current population is used as the reference population. For "VR1", \code{geno} is an allele dosage matrix for bi-allelic loci (sites x indiv). 
#' 
#' For \code{method}="AM", the Allele Matching coefficient is calculated, which is the probability that two haplotypes sampled at random (with replacement) are identical (Weir and Goudet 2017). Missing data are not allowed. For "AM", \code{geno} is a phased genotype matrix, with alleles coded as positive integers. The name of each column is the id and haplotype concatenated with a separating character, \code{sep}. This character needs to be unique (i.e., not present in id or haplotype).
#'
#' @references VanRaden (2008) J. Dairy Sci 91:4414-4423.
#' @references Endelman et al. (2018) Genetics 209:77-87.
#' @references Weir and Goudet (2017) Genetics 206:2085-2103.
#'  
#' @param geno genotype matrix or filename
#' @param ploidy ploidy
#' @param p.ref optional, reference population frequency for method "VR1"
#' @param method "VR1" or "AM"
#' @param sep character separating id from haplotype for "AM" method
#' @param n.core number of cores, only used with "AM"
#' 
#' @return G matrix 
#' 
#' @export
#' @importFrom data.table fread
#' @import parallel

G_mat <- function(geno, ploidy, p.ref=NULL, method="VR1", sep="_", n.core=1) {
  
  stopifnot(method %in% c("VR1","AM"))
  
  if (inherits(geno,"character")) {
    tmp <- fread(geno,header=T,sep=",")
    geno <- as.matrix(tmp[,-1])
    rownames(geno) <- as.character(tmp$V1)
  }
  m <- nrow(geno)
  n <- ncol(geno)
  
  if (method=="VR1") {
    if (is.null(p.ref)) {
      p.ref <- apply(geno,1,mean,na.rm=T)/ploidy
    }
    Pmat <- kronecker(matrix(p.ref,nrow=1,ncol=m),matrix(1,nrow=n,ncol=1))
    coeff <- t(geno) - ploidy * Pmat 
    coeff[is.na(coeff)] <- 0
    dimnames(coeff) <- dimnames(geno)[c(2,1)]
    scale <- ploidy*sum(p.ref*(1-p.ref))
    return(tcrossprod(coeff)/scale)
  } else {
    
    tmp <- strsplit(colnames(geno),split=sep,fixed=T)
    id <- unique(sapply(tmp,"[[",1))
    
    kam <- function(x,ploidy) {
      # x is vector of phased genotypes at one locus
      n <- length(x)/ploidy
      x2 <- matrix(as.character(x),nrow=n,byrow = T)
      haps <- unique(x)
      y <- apply(x2,1,function(u){table(factor(u,levels=haps))})
      if (is.null(dim(y))) {
        y <- matrix(y,nrow=1)
      }
      crossprod(y)/ploidy^2
    }
    geno2 <- split(geno,1:m)
    if (n.core > 1) {
      cl <- makeCluster(n.core)
      clusterExport(cl=cl,varlist="kam",envir=environment())
      u <- parLapply(cl,X=geno2,fun=kam,ploidy=ploidy)
      stopCluster(cl)
    } else {
      u <- lapply(geno2,FUN=kam,ploidy=ploidy)
    }
    G <- ploidy*apply(simplify2array(u),c(1,2),mean)
    dimnames(G) <- list(id,id)
    return(G)
  }
}
