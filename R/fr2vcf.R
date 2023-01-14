#' Convert Genome Studio Final Report to VCF
#' 
#' Convert Genome Studio Final Report to VCF
#'
#' XY values are multiplied by 100 to generate AD field
#' 
#' @param fr.file name of Genome Studio Final Report file
#' @param map.file vcf file with map positions for the markers
#' @param vcf.file output vcf file
#' 
#' @importFrom vcfR read.vcfR
#' @importFrom tidyr pivot_wider
#' 
#' @export

fr2vcf <- function(fr.file, map.file, vcf.file) {
  
  data <- read.table(file=fr.file, skip=9, sep="\t", header=T,as.is=T,check.names=F)
  
  X <- pivot_wider(data=data[,c("SNP Name","Sample ID","X")],
                   names_from="Sample ID",values_from="X")
  id <- colnames(X)[-1]
  marks <- X[,1][[1]]
  Y <- pivot_wider(data=data[,c("SNP Name","Sample ID","Y")],names_from="Sample ID",values_from="Y")
  m <- length(marks)
  n <- length(id)
  tmp <- array(data=0,dim=c(m,n,2))
  tmp[,,1] <- round(as.matrix(X[,-1])*100,0)
  tmp[,,2] <- round(as.matrix(Y[,-1])*100,0)
  AD <- apply(tmp,c(1,2),function(z){
    if(any(is.nan(z)|is.na(z))) {
      return("0,0")
    } else {
      return(paste(z,collapse=","))
    }})
  dimnames(AD) <- list(marks,id)
  
  fAB <- function(x) {
    y <- integer(length(x))
    missing <- (sapply(gregexpr("-",x),function(z){any(z>0)}) | x=="NC")
    if (sum(missing) > 0)
      y[missing] <- NA
    if (sum(!missing) > 0)
      y[!missing] <- sapply(gregexpr("B",x[!missing]),function(z){sum(z>0)})
    return(y)
  }
  
  k <- grep("Allele1 - AB",colnames(data))
  if (length(k)>0) {
    ploidy <- 2
    X <- pivot_wider(data=data[,c("SNP Name","Sample ID","Allele1 - AB")],names_from="Sample ID",values_from="Allele1 - AB")
    Y <- pivot_wider(data=data[,c("SNP Name","Sample ID","Allele2 - AB")],names_from="Sample ID",values_from="Allele2 - AB")
    tmp <- array(data="",dim=c(m,n,2))
    tmp[,,1] <- as.matrix(X[,-1])
    tmp[,,2] <- as.matrix(Y[,-1])
    AB <- apply(tmp,c(1,2),paste,collapse="")
    geno <- apply(AB,2,fAB)
  } else {
    k <- grep("Alleles - AB",colnames(data))
    if (length(k) > 0) {
      AB <- as.matrix(pivot_wider(data=data[,c("SNP Name","Sample ID","Alleles - AB")],names_from="Sample ID",values_from="Alleles - AB")[,-1])
      ploidy <- max(apply(AB,c(1,2),nchar))
      geno <- apply(AB,2,fAB)
    } else {
      geno <- matrix(as.integer(NA),nrow=m,ncol=n)
    }
  }
  dimnames(geno) <- dimnames(AD)
  
  map <- read.vcfR(map.file,verbose=FALSE)
  mark2 <- intersect(marks,map@fix[,"ID"])
  ix <- which(map@fix[,"ID"] %in% mark2)
  mark3 <- map@fix[ix,"ID"]
  
  if (length(mark3)==0) {
    stop("No markers in commmon with  map file")
  }
  GT <- apply(geno[mark3,],c(1,2),make_GT,ploidy=ploidy)
  AD <- AD[mark3,]
  DP <- apply(AD,2,function(x){
      sapply(strsplit(x,split=",",fixed=T),function(u){sum(as.integer(u))})
  })
  
  fixed <- map@fix[ix,]
  fixed[is.na(fixed)] <- "."
  write_vcf(vcf.file, fixed=fixed, geno=list(GT=GT,AD=AD,DP=DP), 
            other.meta=paste0("source=",fr.file))
}
