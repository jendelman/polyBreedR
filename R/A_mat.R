#' Additive relationship matrix from pedigree
#' 
#' Additive relationship matrix from pedigree
#' 
#' This is a wrapper that prepares the pedigree in the format required for R package \code{AGHmatrix} by Amadeu et al. (2016) (cite them if you use this function). A random bivalents model for tetraploid meiosis is assumed.
#' 
#' @references Amadeu et al. (2016) Plant Genome 9, doi:10.3835/plantgenome2016.01.0009
#'  
#' @param ped Pedigree in three column format: id, mother, father
#' @param ploidy 2 or 4
#' @param order.ped TRUE/FALSE does the pedigree need to be ordered so that progeny follow parents
#' 
#' @return Additive relationship matrix (dim: indiv x indiv)
#' 
#' @export
#' @importFrom AGHmatrix Amatrix
#' @importFrom utils capture.output
#' @importFrom pedigree orderPed

A_mat <- function(ped,ploidy,order.ped=TRUE) {
  
  colnames(ped) <- c("id","mother","father")
  for (i in 1:3) {
    ped[,i] <- as.character(ped[,i])
  }
  ped2 <- data.frame(id=1:nrow(ped),
                     parent1=match(ped$mother,ped$id,nomatch=0),
                     parent2=match(ped$father,ped$id,nomatch=0))
  if (order.ped) {
    ix <- orderPed(ped=ped2)
    ped2 <- ped2[order(ix),]
  }
    
  invisible(capture.output(A <- Amatrix(ped2,ploidy=ploidy)))
  rownames(A) <- ped$id
  colnames(A) <- ped$id
  return(A)
}