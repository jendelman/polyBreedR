#' Convert DArTag to VCF
#' 
#' Convert DArTag to VCF
#' 
#' Two input files expected. counts.file is the two-row collapsed counts file, whereas dosage.file has one row per target, with chrom and position in columns 4 and 5. DArT reports dosage of REF, whereas VCF standard is based on dosage of ALT. The dosage is exported as GT field in VCF.
#' 
#' Duplicate samples are renamed by appending the "Target ID".
#' 
#' @param counts.file DArTag collapsed counts file
#' @param dosage.file DArTag dosage file
#' @param vcf.file name of VCF output file (uncompressed)
#' @param ploidy ploidy
#' @param first.data.row default is 9 for DArTag format
#' 
#' @export
#' @importFrom utils read.csv

dart2vcf <- function(counts.file, dosage.file, vcf.file, ploidy,
                     first.data.row=9) {

  data <- read.csv(dosage.file, header=F, check.names=F)
  cols <- 6:ncol(data)
  rows <- first.data.row:nrow(data)
  id <- as.character(data[first.data.row-2,cols])
  dupes <- unique(id[which(duplicated(id))])
  if (length(dupes)>0) {
    cat(paste(c("Duplicates detected for",dupes),collapse="\n"))
    ix <- which(id %in% dupes)
    id[ix] <- apply(cbind(id[ix],as.character(data[first.data.row-1,cols])[ix]),
                    1,paste,collapse=".")
  } 
  
  n <- length(id)
  map <- data[rows,c(1,4,5)]
  m <- nrow(map)
  colnames(map) <- c("marker","chrom","position")
  map$position <- as.integer(map$position)
  geno <- ploidy - apply(data[rows,cols],2,as.integer)
  dimnames(geno) <- list(map$marker,id)
  
  GT <- apply(geno,c(1,2),make_GT,ploidy=ploidy)
  
  data2 <- read.csv(counts.file,header=F,check.names=F)
  cols <- 6:ncol(data2)
  rows <- first.data.row:nrow(data2)
  id <- as.character(data2[first.data.row-2,cols])
  if (length(dupes)>0) {
    ix <- which(id %in% dupes)
    id[ix] <- apply(cbind(id[ix],as.character(data2[first.data.row-1,cols])[ix]),
                    1,paste,collapse=".")
  } 
  
  allele.seq <- matrix(data2[rows,3],ncol=2,byrow = T)
  REF.ALT <- t(apply(allele.seq,1,function(x){
    bases <- c("A","C","G","T")
    hap.ref <- unlist(strsplit(x[1], split = ""))
    hap.alt <- unlist(strsplit(x[2], split = ""))
    k <- which(hap.ref!=hap.alt & hap.ref %in% bases & hap.alt %in% bases)
    if (length(k)==1) {
      return(c(hap.ref[k],hap.alt[k]))
    } else {
      return(c("N","N"))
    }
  }))
  
  counts <- as.matrix(data2[rows,cols])
  
  AD <- apply(counts,2,function(x){
    u <- split(x,rep(1:m,each=2))
    sapply(u,function(v){paste(v,collapse=",")})
  })
  
  marks <- unique(data2[rows,2])
  dimnames(AD) <- list(marks,id)
  rownames(REF.ALT) <- marks
  colnames(REF.ALT) <- c("REF","ALT")

  map <- map[order(map$chrom,map$position),]
  contigs <- unique(map$chrom)
  nchr <- length(contigs)
  
  GT <- GT[map$marker,id,drop=FALSE]
  AD <- AD[map$marker,id,drop=FALSE]
  REF.ALT <- REF.ALT[map$marker,]

  DP <- apply(AD,2,function(x){
    sapply(strsplit(x,split=",",fixed=T),function(u){sum(as.integer(u))})
  })
  
  #put all together
  NS <- apply(DP,1,function(x){sum(x>0)})
  DP.AVG <- apply(DP,1,mean)
  
  alt <- apply(AD,1,function(x){
    sum(as.integer(sapply(strsplit(x,split=",",fixed=T),"[[",2)))
  })
  AF <- round(alt/apply(DP,1,sum),3)
  AF[is.na(AF)] <- "."
  
  info <- make_info(cbind(NS=NS,
                          DP.AVG=round(DP.AVG,1),
                          AF=AF))
  fixed <- cbind(map[,c("chrom","position","marker")],REF.ALT,
        rep(".",m),rep(".",m),info)
  colnames(fixed) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  write_vcf(filename=vcf.file,
            fixed=fixed,
            geno=list(GT=GT,AD=AD,DP=DP))
}


