#' SNP array to VCF
#' 
#' Converts output from Genome Studio (Final Report or Wide) to VCF
#'
#' Auto-detects whether the input file is a Genome Studio Final Report, which is a "long" format with 9-row header, or in "wide" format, where all the data for each marker is one row. XY values are multiplied by 100% to generate AI field. 
#' 
#' Genotype calls will attempt to be imported from the GS Final Report when \code{model.file=NULL}. For diploids, columns named "Allele 1 - AB" and "Allele 2 - AB" are expected. For tetraploids, a single column named "Alleles - AB" is expected. 
#' 
#' It is assumed that the parameters in \code{model.file} lead to genotype calls for the dosage of allele B. For a VCF file, genotype calls need to be based on the dosage of ALT. By default, it is assumed that A is the REF allele. For variants where B is REF, include "REF=B" as INFO in the VCF \code{map.file}. 
#' 
#' @param array.file name of input file with SNP array allele intensities
#' @param map.file vcf file with map positions for the markers
#' @param model.file normal mixture model parameters for genotype calls
#' @param ploidy sample ploidy, for use with \code{model.file}
#' @param vcf.file output vcf file
#' 
#' @importFrom vcfR read.vcfR
#' @importFrom tidyr pivot_wider
#' 
#' @export

array2vcf <- function(array.file, map.file, model.file=NULL, ploidy, vcf.file) {
  
  nc <- nchar(array.file)
  if (substr(array.file,nc-1,nc)==".gz") {
    con <- gzfile(array.file,open="r")
  } else {
    con <- file(array.file,open="r")
  }
  
  tmp <- readLines(con,n=1)
  close(con)
  if (tmp=="[Header]") {
    file.format <- "long"
    data <- read.table(file=array.file, skip=9, sep="\t", header=T,as.is=T,check.names=F)
  
    X <- pivot_wider(data=data[,c("SNP Name","Sample ID","X")],
                     names_from="Sample ID",values_from="X")
    id <- colnames(X)[-1]
    marks <- X$`SNP Name`
    X <- as.matrix(X[,-1])
    
    Y <- pivot_wider(data=data[,c("SNP Name","Sample ID","Y")],
                     names_from="Sample ID",values_from="Y")
    Y <- as.matrix(Y[,-1])
  } else {
    file.format <- "wide"
    data <- read.table(file=array.file, sep="\t", header=T,as.is=T,check.names=F)
    
    tmp <- colnames(data)[-(1:5)]
    tmp <- sub(".R","",tmp,fixed=T)
    tmp <- sub(".Theta","",tmp,fixed=T)
    tmp <- sub(".X","",tmp,fixed=T)
    id <- unique(sub(".Y","",tmp,fixed=T))
    iv <- which(match(paste0(id,".X"),colnames(data),nomatch=0)==0)
    if (length(iv) > 0)
      id <- id[-iv]
    
    X <- as.matrix(data[,paste0(id,".X")])
    Y <- as.matrix(data[,paste0(id,".Y")])
    marks <- data$Name
  }

  m <- length(marks)
  n <- length(id)
  dimnames(X) <- dimnames(Y) <- list(marks,id)
  
  tmp <- array(data=0,dim=c(m,n,2))
  tmp[,,1] <- round(X*100,0)
  tmp[,,2] <- round(Y*100,0)
  
  AI <- apply(tmp,c(1,2),function(z){
    if(any(is.nan(z)|is.na(z))) {
      return("0,0")
    } else {
      return(paste(z,collapse=","))
    }})
  dimnames(AI) <- list(marks,id)
  
  fAB <- function(x) {
    y <- integer(length(x))
    missing <- (sapply(gregexpr("-",x),function(z){any(z>0)}) | x=="NC")
    if (sum(missing) > 0)
      y[missing] <- NA
    if (sum(!missing) > 0)
      y[!missing] <- sapply(gregexpr("B",x[!missing]),function(z){sum(z>0)})
    return(y)
  }
  
  geno <- matrix(as.integer(NA),nrow=m,ncol=n)
  if (is.null(model.file)) {
    if (file.format=="long") {
      #attempt to use genotype calls from Genome Studio Final Report
      k <- grep("Allele1 - AB",colnames(data))
      if (length(k) > 0 & ploidy==2) {
        X <- pivot_wider(data=data[,c("SNP Name","Sample ID","Allele1 - AB")],
                         names_from="Sample ID",values_from="Allele1 - AB")
        Y <- pivot_wider(data=data[,c("SNP Name","Sample ID","Allele2 - AB")],
                         names_from="Sample ID",values_from="Allele2 - AB")
        tmp <- array(data="",dim=c(m,n,2))
        tmp[,,1] <- as.matrix(X[,-1])
        tmp[,,2] <- as.matrix(Y[,-1])
        AB <- apply(tmp,c(1,2),paste,collapse="")
        geno <- apply(AB,2,fAB)
      } 
      
      k <- grep("Alleles - AB",colnames(data))
      if (length(k) > 0 & ploidy==4) {
        #tetraploid
        AB <- as.matrix(pivot_wider(data=data[,c("SNP Name","Sample ID","Alleles - AB")],
                                      names_from="Sample ID",values_from="Alleles - AB")[,-1])
        geno <- apply(AB,2,fAB)
      }
    }
    
    dimnames(geno) <- list(marks,id)
  } else {
    #model.file
    dimnames(geno) <- list(marks,id)
    model.params <- as.matrix(read.csv(file=model.file,as.is=T,row.names = 1))
    mark2 <- intersect(marks,rownames(model.params))
    model.ploidy <- ncol(model.params)/3 - 1
    ratio <- Y/(X+Y)
    dimnames(ratio) <- dimnames(AI)
    geno[mark2,] <- geno_call(data=ratio[mark2,],
                      filename=model.file,
                      model.ploidy=model.ploidy,sample.ploidy=ploidy)
  }
  
  map <- read.vcfR(map.file,verbose=FALSE)
  ix <- which(map@fix[,"ID"] %in% intersect(marks,map@fix[,"ID"]))
  mark3 <- map@fix[ix,"ID"]
  
  if (length(mark3)==0) {
    stop("No markers in commmon with  map file")
  }
  geno <- geno[mark3,]
  AI <- AI[mark3,]
  
  #check for REF information in map file
  iv <- grep("REF=B",map@fix[ix,"INFO"],fixed=T)
  
  if (length(iv) > 0 & !is.null(model.file)) {
    #switch REF and ALT alleles when REF=B
    geno[iv,] <- t(apply(geno[iv,,drop=FALSE],1,function(z){ploidy - z}))
    AI[iv,] <- t(apply(AI[iv,,drop=FALSE],1,function(z){
      tmp <- strsplit(z,split=",",fixed=T)
      sapply(tmp,function(q){paste(q[2],q[1],sep=",")})
    }))
  }
  
  GT <- apply(geno,c(1,2),make_GT,ploidy=ploidy)

  fixed <- map@fix[ix,]
  fixed[is.na(fixed)] <- "."
  
  tmp <- apply(geno,1,function(x){sum(!is.na(x))})*ploidy
  tmp2 <- paste0("AN=",tmp)
  fixed[,"INFO"] <- paste(fixed[,"INFO"],tmp2,sep=";")
  
  tmp <- apply(geno,1,mean,na.rm=T)/ploidy
  tmp2 <- paste0("AF.GT=",ifelse(is.nan(tmp),".",round(tmp,3)))
  fixed[,"INFO"] <- paste(fixed[,"INFO"],tmp2,sep=";")
  
  write_vcf(vcf.file, fixed=fixed, geno=list(GT=GT,AI=AI), 
            other.meta=paste0("source=",array.file))
}
