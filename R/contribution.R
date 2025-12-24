#' Genetic contribution
#' 
#' Computes genetic contribution from pedigree
#'  
#' Argument \code{ped} is three column pedigree; can be created using \code{\link{get_pedigree}}. 
#'  
#' @param id vector of target names
#' @param ped pedigree
#' 
#' @return matrix of contributions, with rows \code{id} and columns for all members of \code{ped}
#' 
#' @export

contribution <- function(id, ped) {

  ped <- as.matrix(ped)
  colnames(ped) <- c("id","mother","father")
  n <- length(id)
  m <- nrow(ped)
  out <- matrix(0,nrow=n,ncol=m,dimnames=list(id,ped[,"id"]))
  i=1
  for (i in 1:n) {
    ix <- match(id[i],ped[,"id"],nomatch = 0)
    if (ix > 0) {
      x <- as.character(na.omit(ped[ix,2:3]))
      gen <- 1
      z <- table(x)*(1/2^gen)
      done <- FALSE
      while(!done) {
        out[i,names(z)] <- out[i,names(z)] + as.numeric(z) #add contribution
        
        #go back another generation
        gen <- gen + 1
        ix <- match(x,ped[,"id"],nomatch = 0)
        if (any(ix > 0)) {
          x <- as.character(na.omit(as.character(ped[ix[ix>0],2:3])))
          z <- table(x)*(1/2^gen)
        } else {
          done <- TRUE
        }
      }
    }
  }
  return(out)
}
