#' Read SNP array intensity data
#' 
#' Read SNP array intensity data
#' 
#' The first two columns of the tab-delimited input file should be the SNP and Sample ID. Columns labeled "X" and "Y" contain the signal intensities for the two alleles. Use \code{output} to specify whether to return the ratio = Y/(X+Y) or theta = atan(Y/X)*2/pi. Option "AD" exports the XY data in the allele depth format for a VCF file ("X,Y"), with the X and Y values multiplied by 100 and rounded to the nearest integer.
#' 
#' @param filename filename 
#' @param skip number of lines to skip before the header line with the column names
#' @param output One of three options: "ratio","theta","AD"
#' 
#' @return matrix with dimensions markers x individuals
#' 
#' @export
#' @importFrom utils read.table
#' @importFrom tidyr pivot_wider
#' 

readXY <- function(filename,skip=9,output="ratio") {
  stopifnot(output %in% c("theta","ratio","AD"))
  
  data <- read.table(file=filename, skip=skip, sep="\t", header=T,as.is=T,check.names=F) 
  X <- pivot_wider(data=data[,c("SNP Name","Sample ID","X")],
                   names_from="Sample ID",values_from="X")
  sample.id <- colnames(X)[-1]
  snp.id <- X[,1][[1]]
  X <- as.matrix(X[,-1])
  Y <- pivot_wider(data=data[,c("SNP Name","Sample ID","Y")],names_from="Sample ID",values_from="Y")
  Y <- as.matrix(Y[,-1])
  
  ADf <- function(z,n) {
    x <- z[1:n]
    y <- z[n+1:n]
    paste(round(100*x,0),round(100*y,0),sep=",")
  }
  data <- switch(output,
                'ratio'=Y/(X+Y),
                'theta'=atan(Y/X)*2/3.14159,
                'AD'=t(apply(cbind(X,Y),1,ADf,n=ncol(X))))
  
  rownames(data) <- snp.id
  colnames(data) <- sample.id
  return(data)
}
