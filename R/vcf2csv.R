#' Convert VCF to CSV
#' 
#' Convert VCF to CSV
#' 
#' @param vcf.file Input file
#' @param csv.file Output file
#' @param format Name of FORMAT key to export, either "GT" or "DS"
#'
#' @return none
#' @export
#' @import vcfR

vcf2csv <- function(vcf.file, csv.file, format) {
  
  stopifnot(format %in% c("GT","DS"))
  data <- read.vcfR(vcf.file,verbose = FALSE)
  meta <- queryMETA(data)
  if (format=="GT") {
    GT <- extract.gt(data,element="GT")
    geno <- apply(GT,2,GT2DS)
  } else {
    stopifnot("FORMAT=ID=DS" %in% meta)
    geno <- extract.gt(data,element="DS",as.numeric=T)
  }
  markers <- data@fix[,"ID"]
  chrom_pos <- apply(data@fix[,c("CHROM","POS")],1,paste,sep="_")
  iv <- is.na(markers)
  markers[iv] <- chrom_pos[iv]
  write.csv(data.frame(marker=markers,chrom=data@fix[,"CHROM"],position=data@fix[,"POS"],geno,check.names = F),file=csv.file,row.names=F)
  
}