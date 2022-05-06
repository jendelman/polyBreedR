#' Process DArTag data
#' 
#' Process DArTag data
#' 
#' Designed for standard two-row format from DArT. Column 1 contains the AlleleID in format MarkerName|Haplotype. Haplotypes are named Ref,RefMatch,Alt,AltMatch,Other. Counts are combined for Ref + RefMatch, as well as Alt + AltMatch. Other haplotypes are discarded. Genotype calls are made using the updog package (Gerard et al. 2018). If a sample has no reads at a marker, the genotype is NA and posterior probability equals 0.
#' 
#' @param filename input filename
#' @param first.data.row first data row
#' @param first.data.col first data column
#' @param ploidy ploidy
#' @param geno.call TRUE/FALSE
#' @param n.core number of cores
#' 
#' @return Data frame with marker statistics and 5 matrices with dimensions markers x indiv
#' \describe{
#' \item{stats}{data frame with 10th percentile of depth and posterior prob for each marker}
#' \item{depth}{(Alt+Ref) read count}
#' \item{ratio}{Alt/(Alt+Ref)}
#' \item{geno.mode}{Posterior mode for Alt allele dosage}
#' \item{prob}{Maximum posterior probability}
#' \item{geno.mean}{Posterior mean for Alt allele dosage}
#' }
#' 
#' @export
#' @importFrom utils read.csv
#' @importFrom updog flexdog
#' @importFrom parallel makeCluster clusterExport parLapply stopCluster
#' @importFrom stats quantile

dart_tag <- function(filename,first.data.row=9,first.data.col=6,ploidy,
                     geno.call=TRUE,n.core=1) {
  
  data <- read.csv(filename,header=F,check.names=F)
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
  
  f1 <- function(x,n,ploidy) {
    tmp <- try(flexdog(refvec=x[1:n],sizevec=x[n+1:n],ploidy=ploidy,model="norm",
                   verbose=FALSE),silent=TRUE)
    if (inherits(tmp,"try-error")) {
      return(c(NA*numeric(2*n),rep(0,n)))
    } 
    ix <- which(x[n+1:n]==0)
    if (length(ix) > 0) {
      tmp$geno[ix] <- NA
      tmp$postmean[ix] <- NA
      tmp$maxpostprob[ix] <- 0
    }
    return(c(tmp$geno,tmp$postmean,tmp$maxpostprob))
  }
  
  out <- list(depth=ref+alt,ratio=alt/(ref+alt))
              
  if (geno.call) {
    tmp <- split(cbind(alt,ref+alt),f=1:m)
  
    if (n.core==1) {
      ans <- lapply(tmp,f1,n=n,ploidy=ploidy)
    } else {
      cl <- makeCluster(n.core)
      clusterExport(cl=cl,varlist=NULL)
      ans <- parLapply(cl=cl,X=tmp,f1,n=n,ploidy=ploidy)
      stopCluster(cl)
    }
    ans <- matrix(unlist(ans),nrow=m,ncol=3*n,byrow = T)
    dimnames(ans) <- list(rownames(alt),rep(colnames(alt),times=3))
    out <- c(out,list(geno.mode=ans[,1:n],
                      prob=ans[,2*n+1:n],geno.mean=ans[,n+1:n]))
    stats <- data.frame(marker=rownames(alt),
                        depth0.1=apply(out$depth,1,quantile,p=0.1),
                        prob0.1=apply(out$prob,1,quantile,p=0.1))
    rownames(stats) <- NULL
  } else {
    stats <- data.frame(marker=rownames(alt),
                        depth0.1=apply(out$depth,1,quantile,p=0.1))
    rownames(stats) <- NULL
  }
  return(c(list(stats=stats),out))
}


