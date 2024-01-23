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
#' @param madc.file MADC filename
#' @param marker Name of marker ("CDF1","OFP20")
#'
#' @return matrix of haplotype counts
#' @export
#' @importFrom data.table fread

madc <- function(madc.file, marker) {

  data <- fread(madc.file,skip=7)
  tmp <- as.matrix(data[,17:ncol(data)])
  id <- colnames(tmp)
  n <- length(id)
  
  if (toupper(marker)=="CDF1") {
    counts <- matrix(0,ncol=4,nrow=n,
                dimnames=list(id,c("CDF1.1","CDF1.2T","CDF1.2C","CDF1.4")))
    k <- which(data$AlleleID=="CDF1.4_chr05_4488021|Ref")
    counts[,"CDF1.1"] <- tmp[k,]
    k <- which(data$AlleleID=="CDF1.4_chr05_4488021|Alt")
    counts[,"CDF1.4"] <- tmp[k,]
    k <- which(data$AlleleID=="CDF1.4_chr05_4488021|Other")
    x <- substr(data$AlleleSequence[k],39,45)
    k2T <- match("ACTAGTG",x,nomatch = 0)
    if (k2T > 0) {
      counts[,"CDF1.2T"] <- tmp[k[k2T],]
    }
    k2C <- match("GCTAGTG",x,nomatch = 0)
    if (k2C > 0) {
      counts[,"CDF1.2C"] <- tmp[k[k2C],]
    }
    k1 <- k[setdiff(1:length(x),c(k2T,k2C))]
    if (length(k1) > 0) {
      counts[,"CDF1.1"] <- counts[,"CDF1.1"] + apply(tmp[k1,,drop=FALSE],2,sum)
    }
    # GT <-apply(counts,1,function(u) {
    #   v <- ifelse(u >= min.AD,TRUE,FALSE)
    #   if (!any(v)) {
    #     return ("3/3")
    #   } else {
    #     w <- c("1","2T","2C","4")[v]
    #     if (length(w) < 2) {
    #       w <- sort(c(w,"3"))
    #     }
    #     paste(w,collapse="/")
    #   }
    # })
    
    return(counts)
  }
  
  if (toupper(marker)=="OFP20") {
    x <- c("OFP20_M6_CDS_994|Ref","OFP20_M6_CDS_994|Alt","OFP20_M6_CDS_24|Ref","OFP20_M6_CDS_24|Alt","OFP20_M6_CDS_171|RefMatch","OFP20_M6_CDS_171|AltMatch")
    ix <- match(x,data$AlleleID)
    data2 <- t(tmp[ix,])
    dimnames(data2) <- list(colnames(tmp),data$AlleleID[ix])
    AF1 <- round(data2[,2]/apply(data2[,1:2],1,sum),2)
    AF8 <- round(data2[,3]/apply(data2[,3:4],1,sum),2)
    data3 <- data.frame(id=rownames(data2),
                        AF1=AF1,
                        AF8=AF8,
                        AD2=data2[,6],
                        'AD3&7'=data2[,5],check.names=F)
                        #AD8=data2[,1])
    rownames(data3) <- NULL
    return(data3)      
  }
}
