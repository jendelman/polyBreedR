#' Generate pedigree
#' 
#' Generate pedigree for a set of individuals
#' 
#' Finds ancestors of individuals in a three-column pedigree file (id,mother,father). The id column can be the identifier for an individual or cross. String matches must be exact or based on the naming convention crossID-progenyID. The returned pedigree is ordered using R package \code{pedigree} so that offspring follow parents. When \code{trim} is TRUE (default), the pedigree is trimmed to remove ancestors with only one offspring (which are not needed to compute the pedigree relationship matrix). 
#'  
#' @param id Vector of names of individuals
#' @param pedfile Name of pedigree file
#' @param delim Delimiter for the pedigree file (default is "," for CSV)
#' @param na.string String used for NA in the pedigree file (default is "NA")
#' @param trim TRUE/FALSE whether to trim pedigree (see Details) 
#' 
#' @return Data frame with columns id, mother, father
#' 
#' @export
#' @importFrom utils read.table
#' @importFrom pedigree orderPed

get_pedigree <- function(id,pedfile,delim=",",na.string="NA",trim=TRUE) {
  
  ped.match <- function(x,y) {
    match2 <- function(x,y) {
      tmp <- sapply(strsplit(x,split="-",fixed=T),function(q){q[1]})
      return(match(tmp,y,nomatch=0))
    }
    return(ifelse(x %in% y,match(x,y,nomatch = 0),match2(x,y)))
  }
  ped <- read.table(file=pedfile, colClasses = rep("character",3), 
                    na.strings=na.string, 
                    sep=delim, header=T)
  colnames(ped) <- c("id","parent1","parent2")
  id.out <- id <- unique(id)
  primary.id <- ped$id
  ped2 <- NULL
    
  while (length(id)>0) {
    ix <- ped.match(id,primary.id)
    present <- which(ix>0)
    missing <- which(ix==0)
    n.miss <- length(missing)
    if (n.miss > 0) {
      tmp <- as.character(rep(NA,n.miss))
      ped2 <- rbind(ped2,data.frame(id=id[missing],parent1=tmp,parent2=tmp,stringsAsFactors=F))
    }
    if (length(present)>0) {
      tmp <- ped[ix[present],]
      ped2 <- rbind(ped2,data.frame(id=id[present],parent1=tmp$parent1,parent2=tmp$parent2,stringsAsFactors = F))
      id <- setdiff(c(tmp$parent1,tmp$parent2),c(ped2$id,NA))  #need to be added
    } else {
      id <- character(0)  #finished
    }
  }

  ped3 <- data.frame(id=1:nrow(ped2),parent1=match(ped2$parent1,ped2$id,nomatch=0),parent2=match(ped2$parent2,ped2$id,nomatch=0))
  rownames(ped3) <- ped2$id
  ix <- orderPed(ped=ped3)
  ped2 <- ped2[order(ix),]
  ped3 <- ped3[order(ix),]
  if (trim) {
    ped3 <- trim_ped(ped3,focal=ped3$id[match(id.out,ped2$id)])
    ped2 <- ped2[ped2$id %in% rownames(ped3),]
  } 
  colnames(ped2) <- c("id","mother","father")
  return(ped2)
}