#' Prepare data for Stage 2 analysis of multi-environment trials
#' 
#' Prepare data for Stage 2 analysis of multi-environment trials
#' 
#' Designed to prepare data files for \code{\link{Stage2}} based on output from \code{\link{Stage1}}. Each element of \code{data} is a list that contains at least two variables: "blue" and "vcov". The "blue" variable is a data frame with columns named "id" and "blue", and if multiple traits have been analyzed in Stage 1, there can be a third column named "trait". The "vcov" variable is the variance-covariance matrix of the BLUEs. By default, the function treats each element of \code{data} as a different environment, which in plant breeding typically represents a location x year combination. If data from multiple experiments per environment are included, each element of \code{data} should also contain the variable "env" to specify the environment name. Furthermore, when the dataset includes multiple locations with more than one environment per location, include "loc" for each element of \code{data} to model genotype x location effects. 
#' 
#' @param data Named list containing output from Stage 1 (see Details)
#' @param id Vector of genotype identifiers to include (default is all)
#' 
#' @return A list containing
#' \describe{
#' \item{blue}{data frame of BLUEs}
#' \item{Omega}{variance-covariance matrix of BLUEs}
#' }
#' For multiple traits, the Omega variable is a list of matrices, one for each trait.
#' 
#' @export
#' @importFrom tidyr pivot_wider

Stage2_prep <- function(data,id=NULL) {
  
  n.expt <- length(data)
  stopifnot(n.expt > 1)
  list.names <- names(data[[1]])
  if (!all(sapply(data,function(x){setequal(names(x),list.names)}))) {
    stop("Each element of data must have the same variables")
  }
  expt <- names(data)
  if (is.null(expt)) {
    stop("data must be a named list")
  }
  
  if ("env" %in% list.names) {
    env <- sapply(data,function(x){x$env})
  } else {
    env <- expt
    expt <- NULL
  }
  if ("loc" %in% list.names) {
    loc <- sapply(data,function(x){x$loc})
    if (min(table(loc)) < 2) {
      stop("Not all locations have at least two environments")
    }
  } else {
    loc <- NULL
  }
  
  data2 <- NULL
  omega <- vector("list",length=n.expt)
  i=1
  for (i in 1:n.expt) {
    tmp <- data.frame(env=env[i],data[[i]]$blue,stringsAsFactors = F)
    tmp$id <- as.character(tmp$id)
    if (!is.null(id)) {
      ix <- which(tmp$id %in% id)
      tmp <- tmp[ix,]
    } else {
      ix <- 1:nrow(tmp)
    }
    if (!is.null(expt)) {tmp$expt <- expt[i]}
    if (!is.null(loc)) {tmp$loc <- loc[i]}
    data2 <- rbind(data2,tmp)
    omega[[i]] <- data[[i]]$vcov[ix,ix]
  }
  Omega <- as.matrix(direct_sum(omega))
  attr(Omega,"INVERSE") <- FALSE
  rownames(Omega) <- colnames(Omega) <- 1:nrow(data2)
  
  if ("trait" %in% colnames(data2)) {
    data2$omega <- 1:nrow(data2)
    data2 <- as.data.frame(pivot_wider(data=data2,names_from="trait",values_from=c("blue","omega")))
    ix <- grep("omega",colnames(data2))
    traits <- sapply(strsplit(colnames(data2)[ix],split="_",fixed=T),function(q){q[2]})
    n.trait <- length(traits)
    omega <- vector("list",n.trait)
    names(omega) <- traits
    for (i in 1:n.trait) {
      iv <- data2[,ix[i]]
      tmp <- Omega[iv,iv]
      rownames(tmp) <- colnames(tmp) <- 1:nrow(tmp)
      attr(tmp,"INVERSE") <- FALSE
      omega[[i]] <- tmp
    }
    data2 <- data2[,-ix]
    ix <- grep("blue",colnames(data2))
    tmp <- colnames(data2)
    tmp[ix] <- traits
    colnames(data2) <- tmp
    Omega <- omega
  }
  return(list(blue=data2,Omega=Omega))
}