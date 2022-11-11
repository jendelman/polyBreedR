#' Convert DArTag file to read count file 
#' 
#' Convert DArTag file to read count file 
#' 
#' Designed for standard two-row format from DArT. Column 1 contains the AlleleID in format MarkerName|Haplotype. Haplotypes are named Ref,RefMatch,Alt,AltMatch,Other. Counts are combined for Ref + RefMatch, as well as Alt + AltMatch. Other haplotypes are discarded. 
#' 
#' Output format contains the read counts for each marker x id combination in the format "count1|count2" (which is the input format for PolyOrigin).
#' 
#' Use \code{AB.file} to convert REF/ALT counts from dart.file to A/B counts from SNP array. The file must have columns named "marker" and "REF", where REF is either A or B. Only markers present in AB.file will be in the output.
#' 
#' @param dart.file DArTag CSV filename
#' @param first.data.row first data row in dart.file
#' @param first.data.col first data column in dart.file
#' @param out.file name of output file
#' @param map.file optional CSV file (marker, chrom, position) to integrate into the output
#' @param AB.file optional CSV file (marker, REF) to convert to allele B dosage
#' 
#' @return marker x indiv matrix of read depths
#' 
#' @export
#' @importFrom utils read.csv

dart_tag <- function(dart.file,first.data.row=9,first.data.col=6,
                     out.file,map.file=NULL,AB.file=NULL) {
  
  data <- read.csv(dart.file,header=F,check.names=F)
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
  
  if (!is.null(AB.file)) {
    AB <- read.csv(AB.file)
    AB <- AB[AB$REF %in% c("A","B"),]
    mark <- intersect(AB$marker,rownames(ref))
    m <- length(mark)
    stopifnot(m > 1)
    ref <- ref[mark,]
    alt <- alt[mark,]
    AB <- AB[match(mark,AB$marker),]
    ib <- AB$marker[AB$REF=="B"]
    tmp <- alt[ib,]
    alt[ib,] <- ref[ib,]
    ref[ib,] <- tmp
  }
  
  ans2 <- apply(rbind(ref,alt),2,
                  FUN=function(u){v <- split(u,f=rep(1:m,times=2))
                  sapply(v,paste,collapse="|")})
  rownames(ans2) <- rownames(alt)
  
  if (!is.null(map.file)) {
    map <- read.csv(map.file)
    colnames(map)[1:3] <- c("marker","chrom","position")
    marks <- intersect(map$marker,rownames(ans2))
    map2 <- map[map$marker %in% marks,1:3]
    write.csv(cbind(map2,ans2[map2$marker,]), file=out.file, row.names=F)
  } else {
    write.csv(ans2,file=out.file)
  }
  return(ref+alt)
}


