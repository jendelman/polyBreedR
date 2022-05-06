#' Convert DArTag data to PolyOrigin format
#' 
#' Convert DArTag data to PolyOrigin format
#' 
#' Designed for standard two-row format from DArT. Column 1 contains the AlleleID in format MarkerName|Haplotype. Haplotypes are named Ref,RefMatch,Alt,AltMatch,Other. Counts are combined for Ref + RefMatch, as well as Alt + AltMatch. Other haplotypes are discarded. Genotype calls can using the updog package (Gerard et al. 2018). If a sample has no reads at a marker, the genotype is NA and posterior probability equals 0.
#' 
#' @param infile input filename
#' @param outfile output filename
#' @param first.data.row first data row
#' @param first.data.col first data column
#' @param array.file optional
#' 
#' @return List of two data frames with statistics for markers and samples
#' 
#' @export
#' @importFrom utils read.csv write.csv
#' @importFrom stats quantile

dart2PO <- function(infile, outfile, first.data.row=9, first.data.col=6, array.file=NULL) {
  
  data <- read.csv(infile,header=F,check.names=F)
  data2 <- as.matrix(data[first.data.row:nrow(data),first.data.col:ncol(data)])
  data2 <- apply(data2,2,as.integer)
  id <- as.character(data[first.data.row-2,first.data.col:ncol(data)])
  
  x <- strsplit(data[first.data.row:nrow(data),1],split="|",fixed=T)
  marker <- sapply(x,"[",1)
  allele <- sapply(x,"[",2)
  allele <- gsub("Match","",allele)
  
  dimnames(data2) <- list(marker,id)

  ans1 <- apply(data2,2,
                FUN=function(u){v <- tapply(u,list(marker,allele),sum)
                                return(c(v[,"Ref"],v[,"Alt"]))})
  m <- nrow(ans1)/2
  n <- ncol(ans1)
  ref <- ans1[1:m,]
  alt <- ans1[m+1:m,]
  if (!is.null(array.file)) {
    map <- read.csv(array.file)
    marks <- intersect(map$marker,rownames(ref))
    map <- map[map$marker %in% marks,]
    B.marks <- map$marker[map$REF=="B"]
    A.marks <- map$marker[map$REF=="A"]
    A <- rbind(ref[A.marks,],alt[B.marks,])
    B <- rbind(ref[B.marks,],alt[A.marks,])
    ref <- A
    alt <- B
  }
  m <- nrow(ref)
  depth <- ref + alt
  
  ans2 <- apply(rbind(ref,alt),2,
                FUN=function(u){v <- split(u,f=rep(1:m,times=2))
                                sapply(v,paste,collapse="|")})
  rownames(ans2) <- rownames(alt)
  write.csv(ans2,file=outfile)
  
  stats1 <- data.frame(marker=rownames(alt),
                      depth0.1=apply(depth,1,quantile,p=0.1))
  rownames(stats1) <- NULL
  
  stats2 <- data.frame(id=colnames(alt),
                       depth0.1=apply(depth,2,quantile,p=0.1))
  rownames(stats2) <- NULL
  
  return(c(list(stats1,stats2)))
}


