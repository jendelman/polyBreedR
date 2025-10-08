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
#' 
#' @return Data frame with columns id, lineage
#' 
#' @export
#' @importFrom utils read.table

maternal_lineage <- function(id, pedfile, delim=",", na.string="NA") {
  
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
  colnames(ped) <- c("id","mother","father")
  n <- length(id)
  ix <- 1:n
  ans <- id
    
  while (length(ix) > 0) {
    iu <- ped.match(ans[ix], ped$id)
    k <- which(iu > 0) 
    if (length(k) > 0) {
      missing <- which(is.na(ped$mother[iu[k]]))
      if (length(missing) > 0)
        k <- k[-missing]
    }
    if (length(k) > 0) {
      ans[ix[k]] <- ped$mother[iu[k]]
      ix <- ix[k]
    } else {
      ix <- integer(0)
    }
  }
  return(data.frame(id=id, lineage=ans))
}