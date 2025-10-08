#' Generate pedigree
#' 
#' Generate pedigree for a set of individuals
#' 
#' Finds ancestors of individuals in a three-column pedigree file (id,mother,father). The id column can be the identifier for an individual or cross. String matches must be exact or based on the naming convention crossID-progenyID. As of v0.48, inbred line names are supported, with en-dash separating each generation of inbreeding. The returned pedigree is ordered using R package \code{pedigree} so that offspring follow parents. When \code{trim} is TRUE (default), the pedigree is trimmed to remove ancestors with only one offspring (which are not needed to compute the pedigree relationship matrix). 
#'  
#' @param id Vector of names of individuals
#' @param pedfile Name of pedigree file
#' @param delim Delimiter for the pedigree file (default is "," for CSV)
#' @param na.string String used for NA in the pedigree file (default is "NA")
#' @param trim TRUE/FALSE whether to trim pedigree (see Details)
#' @param founders TRUE/FALSE should all pedigree founders be included in id column
#' 
#' @return Data frame with columns id, mother, father
#' 
#' @export
#' @importFrom utils read.table
#' @importFrom pedigree orderPed

get_pedigree <- function(id,pedfile,delim=",",na.string="NA",
                         trim=TRUE, founders=F) {
  
  ped.match <- function(x,y) {
    match2 <- function(x,y) {
      tmp <- sapply(strsplit(x,split="-",fixed=T),function(q){q[1]})
      return(match(tmp,y,nomatch=0))
    }
    return(ifelse(x %in% y,match(x,y,nomatch = 0),match2(x,y)))
  }
  
  f.inbred <- function(x) {
    #anything inbred is removed, its parents added
    x2 <- strsplit(x,split="-",fixed=T)
    len <- sapply(x2,length)
    iu <- which(len > 2) 
    n <- length(iu)
    parents <- matrix("",nrow=n,ncol=2)
    if (n > 0) {
      y <- mapply(FUN=function(x,ix){x[ix]},
                    x=x2[iu],ix=as.list(len[iu]-1))
      z <- mapply(FUN=function(x,ix){paste(x[1:ix],collapse="-")},
                    x=x2[iu],ix=as.list(len[iu]-2))
      y <- gsub("(","",y,fixed=T)
      y <- gsub(")","",y,fixed=T)
      for (k in 1:n) {
        sibs <- grep("/",y[k])
        if (length(sibs) > 0) {
          parents[k,] <- paste(z[k],strsplit(y[k],split="/",fixed=T)[[1]],sep="-")
        } else {
          parents[k,] <- paste(z[k],c(y[k],y[k]),sep="-")
        }
      }
    }
    return(list(id=union(setdiff(x,x[iu]),unique(c(parents[,1],parents[,2]))), 
           ped=data.frame(id=x[iu],parent1=parents[,1],parent2=parents[,2])))
  }

  ped <- read.table(file=pedfile, colClasses = rep("character",3), 
                    na.strings=na.string, 
                    sep=delim, header=T)
  colnames(ped) <- c("id","parent1","parent2")
  id.out <- id <- unique(id)
  id <- gsub("-DH","_DH",id)
  primary.id <- ped$id
  ped2 <- NULL
   
  #first deal with inbreds 
  ans <- f.inbred(id)
  while (nrow(ans$ped) > 0) {
    id <- ans$id
    ped2 <- rbind(ped2,ans$ped)
    ans <- f.inbred(id)
  }
  ped2 <- ped2[!duplicated(ped2$id),]
    
  while (length(id)>0) {
    ix <- ped.match(id,primary.id)
    present <- which(ix>0)
    missing <- which(ix==0)
    n.miss <- length(missing)
    if (n.miss > 0) {
      tmp <- as.character(rep(NA,n.miss))
      ped2 <- rbind(ped2,
                    data.frame(id=id[missing],parent1=tmp,parent2=tmp,stringsAsFactors=F))
    }
    if (length(present)>0) {
      tmp <- ped[ix[present],]
      ped2 <- rbind(ped2,
                    data.frame(id=id[present],parent1=tmp$parent1,parent2=tmp$parent2,stringsAsFactors = F))
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
  if (founders) {
    missing <- setdiff(union(ped2$mother,ped2$father),c(ped2$id,as.character(NA)))
    if (length(missing) > 0) {
      ped2 <- rbind(data.frame(id=missing, mother=as.character(NA), father=as.character(NA)), ped2)
    }
  }
  return(ped2)
}