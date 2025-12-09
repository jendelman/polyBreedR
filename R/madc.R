#' Multi-Allelic Haplotype Counts for potato DArTag
#' 
#' Multi-Allelic Haplotype Counts for potato DArTag 
#' 
#' Due to multi-allelism, for some trait markers a correct interpretation is not possible using the collapsed counts file; the MADC (Missing Allele Discovery Count) file is needed.
#' 
#' "CDF1" uses marker CDF1.4_chr05_4488021 to detect the 2C, 2T, and 4 alleles; all other haplotypes are treated as allele 1. Allele 3 is not detected by the assay. 
#' 
#' "OFP20" relies on three markers. Marker OFP20_M6_CDS_994 detects OFP20.1 as Alt and most other haplotypes as Ref, but some alleles appear to be NULL. Marker OFP20_M6_CDS_171 detects allele 2 as Alt and alleles 3 and 7 as Ref; other alleles are NULL. Marker OFP20_M6_CDS_24 detects allele 8 as Ref and most other alleles as Alt. Given the high allelic diversity at this locus, this function may not work in all germplasm groups.
#' 
#' "Rychc" returns read counts at Rychc_M6v5_chr09_37964457 for the M6 allele.
#' 
#' @param madc.file MADC filename
#' @param marker Name of marker ("CDF1","OFP20","Rychc")
#'
#' @return matrix of haplotype counts
#' @export
#' @importFrom data.table fread

madc <- function(madc.file, marker) {

  reverse_complement <- function(x) {
    n <- nchar(x)
    out <- character(n)
    for (i in 1:n) {
      y <- substr(x,i,i)
      k <- match(y,c("A","C","G","T"))
      out[n-i+1] <- switch(k,"T","G","C","A")
    }
    return(paste(out,collapse=""))
  }
  con <- file(madc.file,open="r")
  tmp <- readLines(con,n=10)
  first.row <- grep("AlleleID",tmp)
  ik <- which(strsplit(tmp[1],split=",")[[1]]!="*")
  first.col <- min(ik)
  close(con)
  
  data <- fread(madc.file,skip=first.row-1)
  tmp <- as.matrix(data[,first.col:ncol(data)])
  id <- colnames(tmp)
  n <- length(id)
  
  if (toupper(marker)=="CDF1") {
    counts <- matrix(0,ncol=4,nrow=n,
                dimnames=list(id,c("CDF1.1","CDF1.2T","CDF1.2C","CDF1.4")))
    k <- which(data$CloneID=="CDF1.4_chr05_4488021")
    allele.seq <- apply(array(data$AlleleSequence[k]),1,reverse_complement)
    CDF1.2T <- grep("CAACACTAGTCACTAGG",allele.seq,fixed=T)
    CDF1.2C <- grep("CAACACTAGCCACTAGG",allele.seq,fixed=T)
    CDF1.4 <- grep("CAACACTAGGTATCCCG",allele.seq,fixed=T)
    CDF1.1 <- setdiff(1:length(k),c(CDF1.2T,CDF1.2C,CDF1.4))
    if (!(2 %in% CDF1.4))
      stop("error with CDF1.4")
    if (!(1 %in% CDF1.1))
      stop("error with CDF1.1")
    counts[,"CDF1.1"] <- apply(tmp[k[CDF1.1],,drop=F],2,sum)
    counts[,"CDF1.4"] <- apply(tmp[k[CDF1.4],,drop=F],2,sum)
    if (length(CDF1.2C)>0) 
      counts[,"CDF1.2C"] <- apply(tmp[k[CDF1.2C],,drop=F],2,sum)
    if (length(CDF1.2T)>0) 
      counts[,"CDF1.2T"] <- apply(tmp[k[CDF1.2T],,drop=F],2,sum)

    return(counts)
  }
  
  if (toupper(marker)=="OFP20") {
    x <- c("OFP20_M6_CDS_994|Ref","OFP20_M6_CDS_994|Alt","OFP20_M6_CDS_24|Ref","OFP20_M6_CDS_24|Alt","OFP20_M6_CDS_171|RefMatch","OFP20_M6_CDS_171|AltMatch")
    ix <- match(x,data$AlleleID)
    data2 <- t(tmp[ix,])
    dimnames(data2) <- list(colnames(tmp),data$AlleleID[ix])
    AF1 <- round(data2[,2]/apply(data2[,1:2],1,sum),2)
    AF8 <- round(data2[,3]/apply(data2[,3:4],1,sum),2)
    data3 <- cbind(AF1=AF1,
                    AF8=AF8,
                    AD2=data2[,6],
                   'AD3&7'=data2[,5])
    rownames(data3) <- rownames(data2)
    return(data3)      
  }
  
  if (toupper(marker)=="RYCHC") {
    x <- "Rychc_M6v5_chr09_37964457|Ref"
    ix <- match(x,data$AlleleID)
    data2 <- t(tmp[ix,,drop=FALSE])
    dimnames(data2) <- list(colnames(tmp),"Rychc.M6")
    return(data2)
  }
}
