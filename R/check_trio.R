#' Check markers for parent-offspring trio
#' 
#' Check markers for parent-offspring trio
#' 
#' Computes the percentage of markers at which the two parents and offspring have incompatible allele dosages (for tetraploids, the random bivalents model is used). For dihaploid offspring of a single tetraploid parent, use \code{ploidy} = 4 and "haploid" for the father in \code{parentage}, as well as a diploid (0,1,2) genotype for the offspring. A small amount of missing genotype data can be tolerated.
#' 
#' @param parentage Data frame with three columns: id, mother, father
#' @param geno Matrix of allele dosages: markers x indiv
#' @param ploidy 2 or 4
#' 
#' @return Data frame with the percentage of incompatible markers for each trio
#'
#' @export
#' 
check_trio <- function(parentage,geno,ploidy) {
  stopifnot(ploidy %in% c(2,4))
  if (ploidy==4) {
    zygote <- function(x) {
      if (any(is.na(x))) {
        return(NA)
      } else {
        return(polyBreedR:::zyg4x[x[1]+1,x[2]+1,x[3]+1])
      }
    }
  } else {
    zygote <- function(x) {
      if (any(is.na(x))) {
        return(NA)
      } else {
        return(polyBreedR:::zyg2x[x[1]+1,x[2]+1,x[3]+1])
      }
    }
  }

  colnames(parentage) <- c("id","mother","father")
  for (i in 1:3) {
    parentage[,i] <- as.character(parentage[,i])
  }
  n <- nrow(parentage)
  m <- nrow(geno)
  ans <- numeric(n)*NA
  geno <- cbind(geno,haploid=rep(0,m))

  for (i in 1:n) {
    ix <- match(parentage[i,c(2,3,1)],colnames(geno))
    if (all(!is.na(ix))) {
      ans[i] <- mean(apply(geno[,ix],1,zygote),na.rm=T)
    }
  }
  return(data.frame(parentage,error=round(ans*100,1),stringsAsFactors=F))
}
