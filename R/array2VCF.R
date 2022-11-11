#' Generate pseudo-VCF file from SNP array data
#' 
#' Creates VCF file with GT and AD fields
#' 
#' AD and GT matrices should be markers x indiv. Per the VCF standard, the "." is used for missing data.
#'
#' @param VCF.file name of VCF file to create
#' @param ploidy ploidy
#' @param AD matrix of allele depths, see \code{\link{readXY}}
#' @param GT optional, matrix of genotype dosages (0,1,2..ploidy)
#' @param map optional, 3 column data frame (marker,chrom,pos)
#' @param header optional, header text for the VCF file


array2VCF <- function(VCF.file,ploidy,AD,GT=NULL,map=NULL,header=NULL) {
  
  f1 <- function(x,ploidy) {
    if (is.na(x)) {
      y <- rep(".",ploidy)
    } else {
      stopifnot(x <= ploidy)
      y <- c(rep(0,ploidy-x),rep(1,x))
    }
    return(paste(y,collapse="/"))
  }
  
  out1 <- file(VCF.file,"w")
  id <- colnames(AD)
  n <- length(id)
  m <- nrow(AD)
  markers <- rownames(AD)
  
  if (!is.null(header)) {
    writeLines(con=out1,
               text=paste(c("##",header)))
  }
  
  writeLines(con=out1,
             text=paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",id),
                        collapse="\t"))
  if (!is.null(map)) {
    colnames(map) <- c("marker","chrom","pos")
    chrom <- map$chrom[match(markers,map$marker)]
    chrom[is.na(chrom)] <- "."
    pos <- map$pos[match(markers,map$marker)]
    pos[is.na(pos)] <- "."
  } else {
    chrom <- pos <- rep(".",m)
  }
  
  if (!is.null(GT)) {
    ix <- match(markers,rownames(GT),nomatch = 0)
    iv <- match(colnames(GT),id,nomatch = 0)
  } else {
    ix <- integer(m)
  }
  
  geno1 <- apply(array(rep(as.integer(NA),n)),1,f1,ploidy=ploidy)
  
  for (i in 1:m) {
    geno2 <- geno1
    
    if (ix[i] > 0) {
      geno2[iv] <- apply(array(GT[ix[i],]),1,f1,ploidy=ploidy)
    } 
    tmp <- paste(geno2,AD[i,],sep=":")
    
    writeLines(con=out1,
               text=paste(c(chrom[i],pos[i],markers[i],"A","B",".",".",".","GT:AD",tmp),
                          collapse="\t"))
  }
  close(out1)
  return(0)
}
