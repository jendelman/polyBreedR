#' Create VCFv4.3 file
#' 
#' Create VCFv4.3 file
#' 
#' Several standard INFO key are recognized:
#' ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
#' ##INFO=<ID=AVG.DP,Number=1,Type=Float,Description="Average Sample Depth">
#' ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#' ##INFO=<ID=AB,Number=1,Type=Float,Description="Allelic Bias">
#' ##INFO=<ID=SE,Number=1,Type=Integer,Description="Sequencing Error (PHRED)">
#' ##INFO=<ID=OD,Number=1,Type=Integer,Description="OverDispersion (PHRED)">
#' ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#' ##INFO=<ID=AF.GT,Number=A,Type=Float,Description="Allele Frequency based on GT">
#' ##INFO=<ID=MIN.DP,Number=1,Type=Integer,Description="smallest mean DP for GT group">
#' ##INFO=<ID=HWE.P,Number=1,Type=Integer,Description="p-value for Hardy-Weinberg Equil (PHRED)">
#' ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">"
#' ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">"

#' Every element of \code{geno} is m x n matrix (m variants, n samples), e.g., AD, GT. The FORMAT field is created from the order and names of \code{geno}. Sample names taken from colnames of \code{geno}. Metadata for \code{geno} is generated from the names of the list:
#' ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#' ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele Depth">
#' ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample Depth">
#' ##FORMAT=<ID=DS,Number=1,Type=Float,Description="Posterior Mean Dosage">
#' ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#'
#' Any additional metadata should be included without the ## prefix.
#'
#' @param filename VCF file name
#' @param fixed character matrix with 8 columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
#' @param geno named list of genotype matrices, see Details
#' @param other.meta optional, other metadata (without ##) besides INFO and FORMAT keys
#' 
#' @export

write_vcf <- function(filename, fixed, geno, other.meta=NULL) {

  stopifnot(colnames(fixed)==c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
  con <- file(filename,open="w")
  writeLines(con=con, text="##fileformat=VCFv4.3")
  
  contigs <- apply(array(unique(fixed[,"CHROM"])),1,
              function(z){sub("Q",z,"##contig=<ID=Q>")}) 
  writeLines(con=con, text=contigs)
  
  x <- lapply(strsplit(fixed[,"INFO"],split=";"),strsplit,split="=")
  info.keys <- unique(unlist(lapply(x,function(z){sapply(z,"[[",1)}),recursive=T))
  
  info.meta <- c("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">",
    "##INFO=<ID=AVG.DP,Number=1,Type=Float,Description=\"Average Sample Depth\">",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
    "##INFO=<ID=AB,Number=1,Type=Float,Description=\"Allelic Bias\">",
    "##INFO=<ID=SE,Number=1,Type=Integer,Description=\"Sequencing Error (PHRED)\">",
    "##INFO=<ID=OD,Number=1,Type=Integer,Description=\"OverDispersion (PHRED)\">",
    "##INFO=<ID=AF.GT,Number=1,Type=Float,Description=\"Allele Frequency based on GT\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
    "##INFO=<ID=MIN.DP,Number=1,Type=Integer,Description=\"smallest mean DP for GT group\">",
    "##INFO=<ID=HWE.P,Number=1,Type=Integer,Description=\"p-value for Hardy-Weinberg Equil (PHRED)\">",
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">",
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">")
    
  x2 <- lapply(info.keys,function(x){sub("Q",x,"<ID=Q,")})
  ix <- lapply(x2,grep,x=info.meta,fixed=T)
  iv <- which(sapply(ix,length)==0)
  if (length(iv) > 1)
    stop("Unrecognized keys in INFO")
  
  writeLines(con=con, text=info.meta[unlist(ix)])
  
  m <- nrow(fixed)
  stopifnot(m==nrow(geno[[1]]))
  id <- colnames(geno[[1]])
  n <- length(id)
  if (any(duplicated(id)))
    stop("Duplicate sample names not allowed")
  
  format.keys <- names(geno)
  nf <- length(format.keys)
  if (nf > 1) {
    FORMAT <- paste(format.keys,collapse=":")
  } else {
    FORMAT <- format.keys
  }
  
  format.meta <- c("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">",
  "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample Depth\">",
  "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Posterior Mean Dosage\">",
  "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">")

  ix <- lapply(as.list(paste0("ID=",format.keys)),grep,x=format.meta,fixed=T)
  iv <- which(sapply(ix,length)==0)
  if (length(iv) > 1)
    stop("Unrecognized FORMAT in geno")
  
  writeLines(con=con, text=format.meta[unlist(ix)])
  
  if (!is.null(other.meta))
    writeLines(con=con, text=paste0("##",other.meta))
  
  writeLines(con=con,
             text=paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",id),
                        collapse="\t"))
  
  for (i in 1:m) {
    tmp2 <- character(n)
    for (j in 1:n) {
      tmp <- geno[[1]][i,j]
      for (k in 2:nf) 
        tmp <- c(tmp,geno[[k]][i,j])
      tmp2[j] <- paste(tmp,collapse=":")
    }
    writeLines(con=con,text=paste(c(as.character(fixed[i,]),FORMAT,tmp2),collapse="\t"))
  }
  
  close(con)
}


