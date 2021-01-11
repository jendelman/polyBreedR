#' Stage 1 analysis of multi-environment trials
#' 
#' Stage 1 analysis of multi-environment trials
#' 
#' Stage 1 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation (license is required). The variable \code{data} must have a column labeled "id" with the names of the different genotypes (i.e., clones or individuals). To include other variables (besides "id") in the model, include them in \code{fixed} or \code{random} as appropriate, and make sure they have the correct type in the data frame: factor vs. numeric. If multiple traits are included, a multivariate analysis is performed, and only plots with data for all traits are included. The \code{h2} matrix returned by the function contains the estimated genetic correlations above the diagonal, residual correlations below the diagonal, and plot-based heritability on the diagonal. For multivariate analysis, the data frame \code{blue} returned by the function is in long format, with a column named "trait". By default, the workspace and pworkspace limits for ASReml-R are set at 500mb. If you get an error about insufficient memory, try increasing the appropriate value (workspace for variance estimation and pworkspace for BLUE computation).
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' 
#' @param data Data frame with phenotype data
#' @param traits Vector of column names from \code{data}
#' @param fixed Vector of column names from \code{data} 
#' @param random Vector of column names from \code{data}
#' @param workspace Memory limit for ASRreml-R variance estimation
#' @param pworkspace Memory limit for ASRreml-R BLUE computation
#' @param silent TRUE/FALSE, whether to suppress ASReml-R output
#' 
#' @return List containing
#' \describe{
#' \item{aic}{AIC from ASReml-R}
#' \item{blue}{data frame of BLUEs}
#' \item{vcov}{variance-covariance matrix of the BLUEs}
#' \item{h2}{matrix with heritability, genetic, and residual correlations (see Details)}
#' }
#' 
#' @importFrom stats complete.cases
#' @export

Stage1 <- function(data,traits,fixed=NULL,random=NULL,silent=TRUE,workspace="500mb",pworkspace="500mb") {
  
  stopifnot(requireNamespace("asreml"))
  stopifnot(all(traits %in% colnames(data)))
  if (!is.null(fixed)) {
    stopifnot(all(fixed %in% colnames(data)))
  }
  if (!is.null(random)) {
    stopifnot(all(random %in% colnames(data)))
  }
  if (!is.null(fixed)&!is.null(random)) {
    stopifnot(length(intersect(fixed,random))==0)
  }
    
  data <- data[complete.cases(data[,traits]),]
  data$id <- factor(as.character(data$id))
  
  x <- "asreml(data=data,fixed="
  n.trait <- length(traits)
  if (n.trait > 1) {
    x <- paste0(x,"cbind(",paste(traits,collapse=","),")~")
  } else {
    x <- paste0(x,traits,"~")
  }
  
  n.fixed <- length(fixed)
  if (n.fixed > 0) {
    if (n.trait > 1) {
      fixed2 <- paste(fixed,"trait",sep=":")
    } else {
      fixed2 <- fixed
    }
    x <- paste0(x,paste(fixed2,collapse="+"))
  } else {
    if (n.trait > 1) {
      x <- paste0(x,"trait")
    } else {
      x <- paste0(x,"1")
    }
  }
  
  blup.model <- paste0(x,",")
  if (n.trait > 1) {
    blue.model <- paste0(x,"+id:trait,")
  } else {
    blue.model <- paste0(x,"+id,")
  }
  
  n.random <- length(random)
  if (n.random > 0) {
    if (n.trait > 1) {
      random2 <- paste(random,"idh(trait)",sep=":")
    } else {
      random2 <- random
    }
    x <- paste0("random=~",paste(random2,collapse="+"))
    blue.model <- paste0(blue.model,x,",")
    blup.model <- paste0(blup.model,x,"+")
  } else {
    blup.model <- paste0(blup.model,"random=~")
  }
  if (n.trait > 1) {
    blup.model <- paste0(blup.model,"id(id):us(trait),")
  } else {
    blup.model <- paste0(blup.model,"idv(id),")  
  }
  
  if (n.trait > 1) {
    x <- "residual=~id(units):us(trait))"
  } else {
    x <- "residual=~idv(units))"
  }
  blup.model <- paste0(blup.model,x)
  blue.model <- paste0(blue.model,x)

  cat("Estimating variance components...\n")
  asreml.options(workspace=workspace,pworkspace=pworkspace,maxit=30,trace=!silent)
  blup.ans <- eval(parse(text=blup.model))
  blue.ans <- eval(parse(text=blue.model))
  if (!all(blup.ans$converge,blue.ans$converge)) {
    stop("ASReml-R failed to converge.")
  }
  vc <- summary(blup.ans)$varcomp
  ix <- grep("id",rownames(vc),fixed=T)
  h2 <- matrix(0,nrow=n.trait,ncol=n.trait)
  colnames(h2) <- rownames(h2) <- traits
  if (n.trait==1) {
    iy <- grep("units!units",rownames(vc),fixed=T)
    h2[1,1] <- vc[ix,1]/(vc[ix,1]+vc[iy,1])
  } else {
    Vg <- h2
    Ve <- h2
    tmp <- strsplit(rownames(vc)[ix],split="_",fixed=T)
    tmp <- sapply(tmp,function(z){strsplit(z[2],split=":",fixed=T)[[1]]})
    Vg[cbind(tmp[1,],tmp[2,])] <- vc[ix,1]
    V2 <- Vg
    diag(V2) <- 0
    Vg <- Vg + t(V2)
    h2[upper.tri(h2,diag = F)] <- cov2cor(Vg)[upper.tri(h2,diag = F)]
    
    ix <- grep("units:trait!trait",rownames(vc),fixed=T)
    tmp <- strsplit(rownames(vc)[ix],split="_",fixed=T)
    tmp <- sapply(tmp,function(z){strsplit(z[2],split=":",fixed=T)[[1]]})
    Ve[cbind(tmp[1,],tmp[2,])] <- vc[ix,1]
    V2 <- Ve
    diag(V2) <- 0
    Ve <- Ve + t(V2)
    h2[lower.tri(h2,diag = F)] <- cov2cor(Ve)[lower.tri(h2,diag = F)]
    
    diag(h2) <- diag(Vg/(Vg+Ve))
  }
  
  cat("Computing BLUEs...\n")
  if (n.trait > 1) {
    predans <- predict(blue.ans,classify="id:trait",vcov = TRUE)
    BLUE <- as.data.frame(predans$pvals[,c(1,3,2)])
    colnames(BLUE) <- c("id","blue","trait")
    BLUE$trait <- as.character(BLUE$trait)
  } else {
    predans <- predict(blue.ans,classify="id",vcov = TRUE)
    BLUE <- as.data.frame(predans$pvals[,1:2])
    colnames(BLUE) <- c("id","blue")
  }
  BLUE$id <- as.character(BLUE$id)
  return(list(aic=summary(blue.ans)$aic,blue=BLUE,vcov=predans$vcov,h2=h2))
}
