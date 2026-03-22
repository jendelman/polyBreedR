#' Additive relationship matrix from pedigree
#' 
#' Additive relationship matrix from pedigree
#' 
#' Founders should have NA for mother and father. Any parents not listed in the first column are added as founders. 
#' 
#' This is a wrapper that prepares the pedigree in the format required for R package \code{polyAinv} by Hamilton and Kerr (cite them if you use this function). A random bivalents model for tetraploid meiosis is assumed.
#' 
#' Pedigree assumes mother (dam) before father (sire). 7-col format pedigree is needed for mixed ploidy: id, dam, sire, dam.gamete, dam.lambda, sire.gamete, sire.lambda (Hamilton and Kerr 2018).
#' 
#' @references Hamilton and Kerr (2018) Theor Appl Genet 131:851-860 doi:10.1007/s00122-017-3041-y
#'  
#' @param ped Pedigree in 3 or 7 column format (see Details)
#' @param ploidy integer or NULL for mixed ploidy pedigree
#' @param order.ped TRUE/FALSE does the pedigree need to be ordered so that progeny follow parents
#' 
#' @return Additive relationship matrix (dim: indiv x indiv)
#' 
#' @export
#' @import polyAinv
#' @importFrom utils capture.output
#' @importFrom pedigree orderPed
#' @import Matrix

A_mat <- function(ped, ploidy, order.ped=TRUE) {
  
  if (is.null(ploidy)) {
    stopifnot(ncol(ped)==7)
  } else {
    ploidy <- as.integer(ploidy)
  }
  if (ncol(ped)==3) {
    for (i in 1:3) {
      ped[,i] <- as.character(ped[,i])
    }
    colnames(ped) <- c("id","dam","sire")
    ped <- data.frame(ped, dam.gamete=ploidy/2, dam.lambda=0,
                      sire.gamete=ploidy/2, sire.lambda=0)
  } else {
    colnames(ped) <- c("id","dam","sire","dam.gamete","dam.lambda","sire.gamete","sire.lambda")
  }
  
  parents <- setdiff(union(ped$dam,ped$sire),NA)
  missing.founders <- setdiff(parents,ped$id)
  founders <- ped$id[which(is.na(ped$dam)&is.na(ped$sire))]
  founders.not.parents <- setdiff(founders,parents)
  nfp <- length(founders.not.parents)
  
  if (length(missing.founders) > 0) {
    if (is.null(ploidy))
      stop("Missing founders not allowed for mixed ploidy input")
    ped <- rbind(data.frame(id=missing.founders, 
                            dam=as.character(NA), 
                            sire=as.character(NA),
                            dam.gamete=ploidy/2, dam.lambda=0,
                            sire.gamete=ploidy/2, sire.lambda=0),
               ped)
  }
  if (nfp > 0) {
    ped1 <- ped[-match(founders.not.parents,ped$id),]
  } else {
    ped1 <- ped
  }
  
  ped2 <- data.frame(id=1:nrow(ped1),
                     dam=match(ped1$dam,ped1$id,nomatch=0),
                     sire=match(ped1$sire,ped1$id,nomatch=0),
                     ped1[,4:7])
  if (order.ped) {
    ix <- orderPed(ped=ped2)
    ped2 <- ped2[order(ix),]
  }
    
  # invisible(capture.output(A2 <- Amatrix(ped2,ploidy=ploidy)))
  # rownames(A2) <- ped$id
  # colnames(A2) <- ped$id
  
  invisible(capture.output(ans <- polyAinv(ped2)))
  n <- nrow(ped2)
  Ainv <- sparseMatrix(i=ans$A.inv$id1,
                      j=ans$A.inv$id2,
                      x=ans$A.inv$A.INV,
                      symmetric=TRUE, dims = c(n,n))
  A <- drop0(solve(Ainv),1e-6)
  if (nfp > 0) {
    A <- bdiag(Diagonal(n=nfp),A)
    id <- c(founders.not.parents,ped1$id)
  } else {
    id <- ped1$id
  }
  dimnames(A) <- list(id,id)
  return(as.matrix(A))
}