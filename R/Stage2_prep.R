#' Prepare data for Stage 2 analysis of multi-environment trials
#' 
#' Prepare data for Stage 2 analysis of multi-environment trials
#' 
#' Designed to prepare data files for \code{\link{Stage2}} based on output from \code{\link{Stage1}}. Each element of \code{blue.vcov} is a matrix for one trait, with the first column containing BLUEs and the remainder is their variance-covariance. The rownames are the id and env concatenated with ":". The output \code{blue} is a data frame with id, env, and trait (when multiple traits are provided). This allows for multi-trait analysis in Stage 2. 
#' 
#' @param blue.vcov named list of blue.vcov matrices from Stage 1 
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

Stage2_prep <- function(data,exclude.id=character(0),exclude.env=character(0)) {
  
  traits <- names(data)
  n.trait <- length(data)
  vcov <- vector("list",n.trait)
  blue <- NULL
  for (i in 1:n.trait) {
    tmp <- data[[i]]
    tmp2 <- strsplit(rownames(tmp),split=":",fixed=T)
    id.i <- sapply(tmp2,"[",1)
    env.i <- sapply(tmp2,"[",2)
    ix <- which(!(id.i %in% exclude.id) & !(env.i %in% exclude.env))
    
    blue <- rbind(blue,data.frame(id=id.i[ix],env=env.i[ix],trait=traits[i],blue=as.numeric(tmp[ix,1])))
    vcov[[i]] <- tmp[ix,1+ix]
  }
  
  Omega <- direct_sum(vcov)
  attr(Omega,"INVERSE") <- FALSE
  rownames(Omega) <- colnames(Omega) <- 1:nrow(blue)
  blue$id <- factor(blue$id)
  blue$env <- factor(blue$env)
  blue$trait <- factor(blue$trait)
  if (n.trait==1) {
    blue <- blue[,-match("trait",colnames(blue))] 
  }
  return(list(blue=blue,Omega=Omega))
}