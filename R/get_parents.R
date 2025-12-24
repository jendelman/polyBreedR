#' Get parents of genotypes
#' 
#' Get parents of genotypes
#' 
#' The function generates data frame with the parents of the input \code{id}. For lines with any degree of inbreeding (more than one en-dash) or dihaploids (-DH), no look-up table is needed. Otherwise, the parents are searched in the three-column \code{pedfile} (id,mother,father),  The id column can be the identifier for an individual or cross. String matches must be exact or based on the naming convention crossID-progenyID. 
#'  
#' @param id Vector of target names
#' @param pedfile Name of pedigree file
#' @param delim Delimiter for the pedigree file (default is "," for CSV)
#' @param na.string String used for NA in the pedigree file (default is "NA")
#' @param DH TRUE/FALSE should 4x parent of dihaploid be returned
#' 
#' @return Data frame with columns id, mother, father
#' 
#' @export
#' @importFrom utils read.table

get_parents <- function(id, pedfile, delim=",", na.string="NA", DH=FALSE) {
  
  ped.match <- function(x,y) {
    match2 <- function(x,y) {
      tmp <- sapply(strsplit(x,split="-",fixed=T),function(q){q[1]})
      return(match(tmp,y,nomatch=0))
    }
    return(ifelse(x %in% y,match(x,y,nomatch = 0),match2(x,y)))
  }
  
  ped <- read.table(file=pedfile, colClasses = rep("character",3), 
                    na.strings=na.string, sep=delim, header=T)
  colnames(ped) <- c("id","mother","father")

  n <- length(id)
  output <- data.frame(id=id,
                       mother=as.character(NA),
                       father=as.character(NA))
  
  ix1 <- grep("-DH",id,fixed=T)
  if (length(ix1) > 0) {
    if (DH) {
      output$mother[ix1] <- sapply(strsplit(id[ix1],split="-DH"),"[[",1)
    } #otherwise stays NA
  }
  
  #then inbreds
  ix2 <- setdiff(1:n,ix1)
  if (length(ix2) > 0) {
    x <- strsplit(id[ix2],split="-",fixed=T)
    len <- sapply(x,length)
    iu <- which(len > 2) #some inbreeding
    ix3 <- ix2[which(len <= 2)] #F1
    n2 <- length(iu)
    
    if (n2 > 0) {
      parents <- matrix("",nrow=n2,ncol=2)
      y <- mapply(FUN=function(x,ix){x[ix]},x=x[iu],ix=as.list(len[iu]-1))
      #get the stem
      stem <- mapply(FUN=function(x,ix){paste(x[1:ix],collapse="-")},x=x[iu],ix=as.list(len[iu]-2))
      y <- gsub("(","",y,fixed=T)
      y <- gsub(")","",y,fixed=T)
      for (k in 1:n2) {
        sibs <- grep("/",y[k])
        if (length(sibs) > 0) {
          output[ix2[iu[k]],2:3] <- paste(stem[k],strsplit(y[k],split="/",fixed=T)[[1]],sep="-")
        } else {
          output[ix2[iu[k]],2:3] <- paste(stem[k],c(y[k],y[k]),sep="-")
        }
      }
    }
  }
  
  #now look up remainder in file
  if (length(ix3) > 0) {
    ix <- ped.match(id[ix3],ped$id)
    present <- which(ix>0)
    output[ix3[present],2:3] <- ped[ix[present],2:3]
  }
  
  return(output)
}
