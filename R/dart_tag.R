#' Extract Ref/Alt counts from DArTag data file
#' 
#' Extract Ref/Alt counts from DArTag data file
#' 
#' Designed for standard two-row format from DArT. First 11 rows contain sample information. Column 1 contains the AlleleID in format MarkerName|Haplotype. Haplotypes are named Ref,RefMatch,Alt,AltMatch,Other. Counts are combined for Ref + RefMatch, as well as Alt + AltMatch. Other haplotypes are discarded.
#' 
#' @param filename input filename
#' 
#' @return 3D array of allele counts with dimensions: markers, samples, alleles (ref/alt)
#' @export
#' @importFrom utils read.csv

dart_tag <- function(filename) {
  
  data <- read.csv(filename,header=F,check.names=F)
  x <- strsplit(data[-(1:11),1],split="|",fixed=T)
  marker <- sapply(x,"[",1)
  allele <- sapply(x,"[",2)
  allele <- gsub("Match","",allele)

  v <- as.list(data[-(1:11),-(1:3)])
  ans1 <- apply(data[-(1:11),-(1:3)],2,
                FUN=function(u){v <- tapply(as.integer(u),list(marker,allele),sum)
                                return(c(v[,"Ref"],v[,"Alt"]))})
  m <- nrow(ans1)/2
  n <- ncol(ans1)
  ref <- ans1[1:m,]
  alt <- ans1[m+1:m,]
  id <- as.character(data[7,-(1:3)])
  
  return(array(c(ref,alt),dim=c(m,n,2),
               dimnames=list(rownames(ref),id,c("Ref","Alt"))))
}


