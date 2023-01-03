#' Impute missing marker data 
#' 
#' Impute missing marker data
#' 
#' Assumes input file is sorted by position. Markers with no genetic variance are removed.
#'
#' \code{method="pop"} imputes with the population mean for \code{geno="DS"} and population mode for \code{geno="GT"}.
#' 
#' \code{method="EM"} uses parameter "tol" (default is 0.02, see rrBLUP A.mat documentation). Imputed values are truncated if needed to fall in the interval [0,ploidy].
#'  
#' \code{method="RF"} uses parameters "ntree" (default 100) for number of trees and "nflank" (default 100) for the number of flanking markers (on each side) to use as predictors. Because RF first uses EM to generate a complete dataset, parameter "tol" is also recognized.
#'
#' 
#' @param in.file VCF input file
#' @param out.file VCF output file
#' @param ploidy ploidy
#' @param method One of the following: "pop","EM","RF"
#' @param geno One of the following: "GT","DS"
#' @param min.DP genotypes below this depth are set to missing (default=1)
#' @param max.missing remove markers above this threshold, as proportion of population
#' @param params list of method-specific parameters
#' @param n.core multicore processing 
#' 
#' @export
#' @import vcfR
#' @importFrom rrBLUP A.mat
#' @importFrom parallel makeCluster stopCluster parApply clusterExport
#' @importFrom randomForest randomForest
#' @importFrom stats sd
#' @importFrom utils read.csv write.csv

impute <- function(in.file,out.file,ploidy,method,geno,min.DP=1,
                   max.missing,params=NULL,n.core=1) {
  
  stopifnot(geno %in% c("GT","DS"))
  stopifnot(method %in% c("pop","EM","RF"))
  stopifnot(min.DP >= 1)

  if (!is.null(params) && !is.list(params))
    stop("Use a list for argument params")
  
  if (method=="EM")
    params1 <- list(tol=0.02)
  if (method=="RF")
    params1 <- list(ntree=100,nflank=100,tol=0.02)
  
  if (!is.null(params)) {
    if (method=="EM")
      stopifnot(names(params)=="tol")
    if (method=="RF")
      stopifnot(names(params) %in% c("ntree","nflank","tol"))
    
    params1[names(params)] <- params
  } 
  
  impute.mode <- function(x) {
    miss <- which(is.na(x))
    if (length(miss)>0) {
      x[miss] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  impute.mean <- function(x) {
    miss <- which(is.na(x))
    if (length(miss)>0) {
      x[miss] <- mean(x[-miss])
    }
    return(x)
  }
  convert <- function(x,field) {
    if (field=="DS")
      return(as.numeric(x))
    if (field=="GT") {
      z <- as.integer(strsplit(x,split="/",fixed=T)[[1]])
      return(sum(z))
    }
  }
  impute.RF <- function(y,x,ntree) {
    miss <- which(is.na(y))
    n.miss <- length(miss)
    if (n.miss > 0) {
      ans <- suppressWarnings(randomForest(x=x[-miss,],
                          y=y[-miss],
                          xtest=x[miss,,drop=FALSE],
                          ntree=ntree))
      y[miss] <- ans$test$predicted
    } 
    if (is.factor(y)) {
      return(as.integer(as.character(y)))
    } else {
      return(y)
    }
  }

  data <- read.vcfR(file=in.file)
  fields <- vcf_field_names(data, tag = "FORMAT")$ID
  geno.key <- geno
  if (!(geno.key %in% fields))
    stop(sub("Q",geno.key,"Q not present in VCF file"))
  
  geno <- extract.gt(data,element=geno.key)
  if (min.DP > 1) {
    dnames <- dimnames(geno)
    geno <- split(geno,f=1:nrow(geno))
    if (!("DP" %in% fields))
      stop("Sample DP not present in VCF file")
    DP <- lapply(split(extract.gt(data,element="DP"),f=1:length(geno)),as.integer)
    geno <- t(mapply(FUN=function(x,dp,min.DP){x[dp < min.DP] <- NA; x},geno,DP,min.DP=min.DP))
    dimnames(geno) <- dnames
  }
  geno <- apply(geno,MARGIN=c(1,2),FUN=convert,field=geno.key)
  
  n <- ncol(geno)
  ix <- which(apply(geno,1,sd,na.rm=T) > 0)
  nd <- nrow(geno) - length(ix)
  if (nd > 0)
    cat(sub("X",nd,"Removed X markers without genetic variance\n"))
  
  iu <- which(apply(geno[ix,],1,function(x){sum(is.na(x))})/n <= max.missing)
  nd2 <- length(ix) - length(iu)
  if (nd2 > 0)
    cat(sub("X",nd2,"Removed X markers due to missing data\n"))
  
  geno <- geno[ix[iu],,drop=FALSE]
  m <- nrow(geno)
  marks <- rownames(geno)
  
  #begin imputation
  #geno.imp is transposed vs geno
  if (method=="pop") {
    if (geno.key=="GT") {
      geno.imp <- apply(geno,1,impute.mode)
    } else {
      geno.imp <- apply(geno,1,impute.mean)
    }
  }
  
  if (method %in% c("EM","RF")) {
    ans <- A.mat(t(geno)/(ploidy/2) - 1,
                 impute.method = "EM",n.core=n.core,
                 return.imputed=TRUE,min.MAF=0,tol=params1$tol)
    geno.imp <- apply(ans$imputed,2,function(x){ifelse(abs(x)>1,1*sign(x),x)})
    if (geno.key=="GT") {
      digits <- 0
    } else {
      digits <- 2
    }
    geno.imp <- round((geno.imp+1)*ploidy/2,digits)
  }
  
  if (method=="RF") {
    geno.imp.old <- geno.imp
    if (n.core > 1) {
      cl <- makeCluster(n.core)
      clusterExport(cl=cl,varlist=NULL)
      geno.imp <- parApply(cl,array(1:m),MARGIN=1,
                      function(k){
                        y <- geno[k,]
                        if (geno.key=="GT")
                          y <- factor(y)
                        impute.RF(y=y,
                        x=geno.imp.old[,setdiff(c(max(1,k-params1$nflank):min(m,k+params1$nflank)),k),drop=FALSE],
                        ntree=params1$ntree)})
      stopCluster(cl)
    } else {
      geno.imp <- apply(array(1:m),MARGIN=1,
                      function(k){
                        y <- geno[k,]
                        if (geno.key=="GT")
                          y <- factor(y)
                        impute.RF(y=y,
                        x=geno.imp.old[,setdiff(c(max(1,k-params1$nflank):min(m,k+params1$nflank)),k),drop=FALSE],
                        ntree=params1$ntree)})
    }
  }
  
  if (geno.key=="GT") {
    geno.imp <- apply(t(geno.imp),c(1,2),make_GT,ploidy=ploidy)
  } else {
    geno.imp <- round(t(geno.imp),1)
  }
  dimnames(geno.imp) <- dimnames(geno)
  
  if (geno.key=="GT")
    fields <- setdiff(fields,"GQ")
  nf <- length(fields)
  geno <- vector("list",nf)
  names(geno) <- fields
  for (i in 1:nf) {
    if (fields[i]==geno.key) {
      geno[[i]] <- geno.imp
    } else {
      geno[[i]] <- extract.gt(data,element=fields[i])[marks,]
    }
  }
  fixed <- data@fix[match(marks,data@fix[,"ID"]),]
  fixed[,"INFO"] <- "."
  fixed[is.na(fixed)] <- "."
  write_vcf(filename=out.file, fixed=fixed, geno=geno)
}