#' Multi-Allelic Haplotype Counts from DArTag
#' 
#' Multi-Allelic Haplotype Counts from the DArTag MADC (Missing Allele Discovery Count) file
#' 
#' Due to multi-allelism, for some trait markers a correct interpretation is not possible using the collapsed counts file; the MADC file is needed. Currently, the only marker implemented is CDF1 for potato DArTag. The CDF1 marker detects the 2C, 2T, and 4 alleles, and all other haplotypes are treated as allele 1. Allele 3 is not detected by the assay. 
#' 
#' @param madc.file MADC filename
#' @param marker Name of marker ("CDF1" is only option so far)
#'
#' @return matrix of haplotype counts
#' @export
#' @importFrom data.table fread

madc <- function(madc.file, marker="CDF1") {
  
  if (toupper(marker)=="CDF1") {
    data <- fread(madc.file,skip=7)
    tmp <- as.matrix(data[,17:ncol(data)])
    id <- colnames(tmp)
    n <- length(id)
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
}
