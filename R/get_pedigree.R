#' Generate pedigree 
#' 
#' Generate pedigree for A matrix calculation
#' 
#' Returns ancestors of \code{id} using lower-level function \code{\link{get_parents}}. Input \code{pedfile} is three-column pedigree file (id,mother,father). The id column can be the identifier for an individual or cross. String matches must be exact or based on the naming convention crossID-progenyID. As of v0.48, inbred line names are supported, with en-dash separating each generation of inbreeding. 
#' 
#' The returned pedigree is ordered using R package \code{pedigree} so that offspring follow parents. When \code{trim} is TRUE (default), the pedigree is trimmed to remove ancestors with only one offspring (which are not needed to compute the pedigree relationship matrix). 
#'  
#' @param id Vector of genotype names
#' @param pedfile Name of pedigree file
#' @param delim Delimiter for the pedigree file (default is "," for CSV)
#' @param na.string String used for NA in the pedigree file (default is "NA")
#' @param trim TRUE/FALSE whether to trim pedigree (see Details)
#' @param founders TRUE/FALSE should all pedigree founders be included in id column
#' @param DH TRUE/FALSE should 4x parent of dihaploids be included
#' 
#' @return Data frame with columns id, mother, father
#' 
#' @export
#' @importFrom utils read.table
#' @importFrom pedigree orderPed

get_pedigree <- function(id,pedfile,delim=",",na.string="NA",
                         trim=TRUE, founders=F, DH=F) {
  
  ped <- read.table(file=pedfile, colClasses = rep("character",3), 
                    na.strings=na.string, 
                    sep=delim, header=T)
  colnames(ped) <- c("id","mother","father")
  
  ped2 <- get_parents(id,pedfile,delim=delim,na.string=na.string,DH=DH)
  new.id <- setdiff(union(ped2$mother,ped2$father),
                    c(ped2$id,as.character(NA))) 
  while (length(new.id) > 0) {
    tmp <- get_parents(new.id,pedfile,delim=delim,na.string=na.string,DH=DH)
    ped2 <- rbind(tmp,ped2)
    new.id <- setdiff(union(ped2$mother,ped2$father),
                      c(ped2$id,as.character(NA))) 
  }
  iv <- which(!duplicated(ped2$id))
  ped2 <- ped2[iv,]
  ped3 <- data.frame(id=1:nrow(ped2),
                     mother=match(ped2$mother,ped2$id,nomatch=0),
                     father=match(ped2$father,ped2$id,nomatch=0))
  rownames(ped3) <- ped2$id
  ix <- orderPed(ped=ped3)
  ped2 <- ped2[order(ix),]
  ped3 <- ped3[order(ix),]
  if (trim) {
    ped3 <- trim_ped(ped3,focal=ped3$id[match(id,ped2$id)])
    ped2 <- ped2[ped2$id %in% rownames(ped3),]
  } 
  
  if (founders) {
    missing <- setdiff(union(ped2$mother,ped2$father),
                       c(ped2$id,as.character(NA)))
    if (length(missing) > 0) {
      ped2 <- rbind(data.frame(id=missing, 
                               mother=as.character(NA), 
                               father=as.character(NA)), ped2)
    }
  }
  return(ped2)
}