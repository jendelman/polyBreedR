#' Impute from low to high density markers by Random Forest
#' 
#' Impute from low to high density markers by Random Forest
#' 
#' Argument \code{params} is a list with three named elements: format, n.tree, n.mark. \code{format} can have values "GT" (integer dosage) or "DS" (real numbers between 0 and ploidy). Classification trees are used for GT and regression trees for DS. \code{n.tree} is the number of trees (default = 100). \code{n.mark} is the number of markers to use as predictors (default = 200), chosen based on minimum distance to the target. 
#' 
#' The \code{exclude} argument is useful for cross-validation. 
#' 
#' Both VCF and CSV are allowable input file formats--they are recognized based on the file extension. For CSV, the first three columns should be marker, chrom, pos. The output file is CSV. 
#' 
#' Any missing data are imputed separately for each input file at the outset, using the population mean (DS) or mode (GT) for each marker. 
#' 
#' @param high.file name of high density file
#' @param low.file name of low density file
#' @param out.file name of CSV output file for imputed data
#' @param params list of parameters (see Details)
#' @param exclude optional, vector of high density samples to exclude
#' @param n.core multicore processing 
#' 
#' @return matrix of OOB error with dimensions markers x trees
#' 
#' @import vcfR
#' @importFrom stats cov sd na.omit
#' @importFrom parallel makeCluster stopCluster parApply clusterExport
#' @importFrom randomForest randomForest
#' @importFrom data.table fread fwrite
#' @export

impute_L2H <- function(high.file, low.file, out.file, params=list(), 
                        exclude=NULL, n.core=1) {
 
  impute.mode <- function(x) {
    miss <- which(is.na(x))
    if (length(miss)>0) {
      x[miss] <- names(which.max(table(x)))
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
  impute.RF <- function(data, geno, train.id, pred.id, params) {
    y <- data$y
    ntest <- nrow(geno) - length(train.id)
    if (params$format=="GT") {
      if (length(unique(y))==1) {
        return(list(pred=rep(y[1],ntest), error=rep(0,params$n.tree)))
      } else {
        if (is.null(data$pred)) {
          imp <- as.integer(names(which.max(table(y))))
          return(list(pred=rep(imp,ntest), error=rep(as.numeric(NA),params$n.tree)))
        } else {
          ans <- randomForest(y=factor(y), x=geno[train.id,data$pred,drop=FALSE],
                            xtest=geno[pred.id,data$pred,drop=FALSE],
                            ntree=params$n.tree)
          return(list(pred=as.integer(as.character(ans$test$predicted)), error=ans$err.rate[,1]))
        }
      }
    } else {
      if (sd(y)==0) {
        return(list(pred=rep(y[1],ntest), error=rep(0,params$n.tree)))
      } else {
        if (is.null(data$pred)) {
          imp <- mean(y,na.rm=T)
          return(list(pred=rep(imp,ntest), error=rep(as.numeric(NA),params$n.tree)))
        } else {
          ans <- suppressWarnings(randomForest(y=y, x=geno[train.id,data$pred,drop=FALSE],
                                             xtest=geno[pred.id,data$pred,drop=FALSE],
                                             ntree=params$n.tree))
          return(list(pred=ans$test$predicted, error=ans$mse))
        }
      }
    }
  }
  
  np <- names(params)
  if (!("format" %in% np))
    params$format <- "GT"
  if (!("n.tree" %in% np))
    params$n.tree <- 100
  if (!("n.mark" %in% np))
    params$n.mark <- 200
  
  vcf1 <- length(grep("VCF",toupper(high.file),fixed=T)) > 0
  if (vcf1) {
    high <- read.vcfR(file=high.file,verbose=FALSE)
    
    fields <- vcf_field_names(high, tag = "FORMAT")$ID
    if (!(params$format %in% fields))
      stop("Invalid FORMAT")
  
    if (params$format=="GT") {
      #geno1 is oriented id x marker
      geno1 <- t(GT2DS(apply(extract.gt(high,element="GT"), 2, impute.mode),
                   n.core=n.core))
    } else {
      geno1 <- t(apply(extract.gt(high,element="DS",as.numeric=T), 2, impute.mean))
    }
    
    map1 <- data.frame(marker=high@fix[,"ID"], 
                       chrom=high@fix[,"CHROM"], 
                       pos=as.integer(high@fix[,"POS"]))
    map1$marker[which(is.na(map1$marker) | map1$marker==".")] <- paste(map1$chrom,map1$pos,sep="_")
    
  } else {
    high <- fread(high.file,check.names=F) 
    map1 <- high[,1:3]
    colnames(map1) <- c("marker","chrom","pos")
    geno1 <- t(as.matrix(high[,-(1:3)]))
    if (params$format=="GT")
      geno1 <- apply(geno1,2,as.integer)
    dimnames(geno1) <- list(colnames(high)[-(1:3)],map1$marker)
  }
  if ("exclude" %in% np) {
    geno1 <- geno1[setdiff(rownames(geno1),params$exclude),]
  }
  
  vcf2 <- length(grep("VCF",toupper(low.file),fixed=T)) > 0
  if (vcf2) {
    low <- read.vcfR(file=low.file,verbose=FALSE)
    
    fields <- vcf_field_names(low, tag = "FORMAT")$ID
    if (!(params$format %in% fields))
      stop("Invalid FORMAT")
    
    if (params$format=="GT") {
      geno2 <- t(GT2DS(apply(extract.gt(low,element="GT"), 2, impute.mode),
                       n.core=n.core))
    } else {
      geno2 <- t(apply(extract.gt(low,element="DS",as.numeric=T), 2, impute.mean))
    }
    
    map2 <- data.frame(marker=low@fix[,"ID"], 
                       chrom=low@fix[,"CHROM"], 
                       pos=as.integer(low@fix[,"POS"]))
    map2$marker[which(is.na(map2$marker) | map2$marker==".")] <- paste(map2$chrom,map2$pos,sep="_")
    
  } else {
    low <- fread(low.file,check.names=F) 
    map2 <- low[,1:3]
    colnames(map2) <- c("marker","chrom","pos")
    geno2 <- t(as.matrix(low[,-(1:3)]))
    if (params$format=="GT")
      geno2 <- apply(geno2,2,as.integer)
    dimnames(geno2) <- list(colnames(low)[-(1:3)],map2$marker)
  }
  
  train.id <- intersect(rownames(geno1),rownames(geno2))
  pred.id <- setdiff(rownames(geno2),train.id)
  if (length(pred.id)==0)
    stop("All low-density individuals in high-density file.")
    
  nid <- length(train.id)
  cat(sub("Q",nid,"Training Set N=Q\n"))
  if (nid < 5) {
    stop("Insufficient Training Set")
  }
  geno1 <- geno1[train.id,]
  
  m1 <- ncol(geno1)
  m2 <- ncol(geno2)
  
  if (n.core > 1) {
    cl <- makeCluster(n.core)
    clusterExport(cl=cl,varlist=NULL)
  }
  
  data <- vector("list",m1)
  names(data) <- map1$marker
  for (i in 1:m1) {
    ix <- which(map2$chrom==map1$chrom[i])
    zz <- min(params$n.mark,length(ix))
    if (zz > 0) {
      data[[i]] <- list(y=geno1[,i],
                 pred=map2$marker[ix[order(abs(map2$pos[ix]-map1$pos[i]))[1:zz]]])
    } else {
      data[[i]] <- list(y=geno1[,i], pred=NULL)
    }
  }
  
  if (n.core > 1) {
    tmp <- parLapply(cl, X=data, fun=impute.RF, geno=geno2, train.id, pred.id, params)
  } else {
    tmp <- lapply(data, FUN=impute.RF, geno=geno2, train.id, pred.id, params)
  }
  geno.imp <- matrix(sapply(tmp,"[[",1),nrow=m1,byrow = T)
  dimnames(geno.imp) <- list(map1$marker,pred.id)
  geno.out <- cbind(t(geno1[train.id,]),geno.imp)[,rownames(geno2)]
  error <- matrix(sapply(tmp,"[[",2),nrow=m1,byrow = T)
  
  if (n.core > 1) {
    stopCluster(cl)
  }
  
  # if (vcf1) {
  #   GT <- apply(geno.out,c(1,2),make_GT,ploidy=ploidy)
  #   info <- make_info(cbind(NS=ncol(geno.out),
  #                           OOB=error[,params$n.tree]))
  #   fixed <- cbind(high@fixed[,1:5],
  #                  rep(".",m1),rep(".",m1),info)
  #   colnames(fixed) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  #   write_vcf(filename=out.file,
  #             fixed=fixed,
  #             geno=list(GT=GT,AD=AD,DP=DP))
  # } else {
  fwrite(cbind(map1,geno.out),file=out.file,row.names=F)
  
  return(error)
}
