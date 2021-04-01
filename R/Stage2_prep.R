#' Prepare data for single-trait, Stage 2 analysis of multi-environment trials
#' 
#' Prepare data for single-trait, Stage 2 analysis of multi-environment trials
#' 
#' Designed to prepare data for \code{\link{Stage2}} based on output from \code{\link{Stage1}}. The first column of \code{blue.vcov} contains BLUEs, and the remainder is their variance-covariance. The rownames are the id and env concatenated with ":". 
#' 
#' @param blue.vcov blue.vcov matrix from Stage 1 
#' @param trait trait name (used to name column in output)
#' @param exclude.id vector of individuals to exclude
#' @param exclude.env vector of envs to exclude
#' 
#' @return a list containing
#' \describe{
#' \item{blue}{data frame of BLUEs}
#' \item{Omega}{variance-covariance matrix of BLUEs}
#' }
#' 
#' @export

Stage2_prep <- function(blue.vcov,trait,
                        exclude.id=character(0),exclude.env=character(0)) {
  
  blue.vcov <- unlist(blue.vcov)
  tmp <- strsplit(rownames(blue.vcov),split=":",fixed=T)
  id.i <- sapply(tmp,"[",1)
  env.i <- sapply(tmp,"[",2)
  ix <- which(!(id.i %in% exclude.id) & !(env.i %in% exclude.env))
  blue <- data.frame(id=id.i[ix],env=env.i[ix],blue=as.numeric(blue.vcov[ix,1]))
  colnames(blue) <- c("id","env",trait)
  Omega <- blue.vcov[ix,1+ix]
  attr(Omega,"INVERSE") <- FALSE
  rownames(Omega) <- colnames(Omega) <- 1:nrow(Omega)
  return(list(blue=blue,Omega=Omega))
  
}

