#' Impute from low to high density markers
#' 
#' Impute from low to high density markers
#' 
#' Missing data in high density file is imputed with the population mean (DS) or mode (GT) for each marker. The Random Forest package is used for imputing up to high density. Argument \code{params} recognizes the following parameters: "format","n.tree","n.mark". Format can have values "GT" or "DS" (default = "GT"), n.tree is the number of trees (default = 100), and n.mark is the number of markers to use as predictors (default = 200), chosen based on minimum distance to the target. Classification trees are used for GT and regression trees for DS. Any markers in low.file that are not present in the high.file or with an identical position are discarded.
#' 
#' @param high.file name of high density VCF file
#' @param low.file name of low density VCF file
#' @param out.file name of CSV output file for the imputed dataset
#' @param ploidy ploidy
#' @param params list of parameters (see Details)
#' @param n.core multicore processing 
#' 
#' @return vector of mean OOB error (vs. number of trees)
#' 
#' @import vcfR
#' @importFrom stats cov sd
#' @importFrom parallel makeCluster stopCluster parApply clusterExport
#' @importFrom randomForest randomForest
#' @importFrom utils write.csv
#' @export

impute_L2H <- function(high.file, low.file, out.file, ploidy,
                        params=list(format="GT", n.tree=100, n.mark=100),
                        n.core=1) {
  
  stopifnot(c("format","n.tree","n.mark") %in% names(params))
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
  impute.RF <- function(i,geno.train,geno.test,pred.ix,params) {
    y <- geno.train[,i]
    ntest <- nrow(geno.test)
    if (params$format=="GT") {
      if (length(unique(y))==1) {
        return(list(pred=rep(y[1],ntest), error=rep(0,params$n.tree)))
      } else {
        ans <- randomForest(y=factor(y), x=geno.train[,pred.ix[i,]],
                            xtest=geno.test[,pred.ix[i,]],
                            ntree=params$n.tree)
        return(list(pred=as.character(ans$test$predicted), error=ans$err.rate[,1]))
      }
    } else {
      if (sd(y)==0) {
        return(list(pred=rep(y[1],ntest), error=rep(0,params$n.tree)))
      } else {
        ans <- suppressWarnings(randomForest(y=y, x=geno.train[,pred.ix[i,]],
                                             xtest=geno.test[,pred.ix[i,]],
                                             ntree=params$n.tree))
        return(list(pred=ans$test$predicted, error=ans$mse))
      }
    }
  }
  
  data <- read.vcfR(file=high.file,verbose=FALSE)
  map1 <- data.frame(marker=data@fix[,"ID"], 
                     chrom=data@fix[,"CHROM"], 
                     pos=as.integer(data@fix[,"POS"]))
  map1$marker[which(is.na(map1$marker) | map1$marker==".")] <- paste(map1$chrom,map1$pos,sep=":")
  m1 <- nrow(map1)
  cp1 <- paste(map1$chrom, map1$pos, sep=":")
  
  fields <- vcf_field_names(data, tag = "FORMAT")$ID
  if (!(params$format %in% fields))
    stop("Invalid FORMAT")
  
  if (params$format=="GT") {
    geno1 <- apply(apply(extract.gt(data,element="GT"), 2, impute.mode),1,GT2DS) #transposes
    rownames(geno1) <- colnames(data@gt)[-1]
  } else {
    geno1 <- t(apply(extract.gt(data,element="DS",as.numeric=T), 2, impute.mean))
  }
  colnames(geno1) <- data@fix[,"ID"]

  data <- read.vcfR(file=low.file,verbose=FALSE)
  map2 <- data.frame(marker=data@fix[,"ID"], 
                     chrom=data@fix[,"CHROM"], 
                     pos=as.integer(data@fix[,"POS"]))
  map2$marker[which(is.na(map2$marker) | map2$marker==".")] <- paste(map2$chrom,map2$pos,sep=":")
  
  fields <- vcf_field_names(data, tag = "FORMAT")$ID
  if (!(params$format %in% fields))
    stop("Invalid FORMAT")
  
  if (params$format=="GT") {
    geno2 <- apply(apply(extract.gt(data,element="GT"), 2, impute.mode),1,GT2DS)
    rownames(geno2) <- colnames(data@gt)[-1]
  } else {
    geno2 <- t(apply(extract.gt(data,element="DS",as.numeric=T), 2, impute.mean))
  }
  colnames(geno2) <- data@fix[,"ID"]
  
  map2 <- map2[map2$marker %in% map1$marker,]
  m2 <- nrow(map2) 
  if (m2==0)
    stop("There are no common markers.")
  geno2 <- geno2[,map2$marker]
  
  map.imp <- map1[!(map1$marker %in% map2$marker),]
  mi <- nrow(map.imp)
  if (mi==0)
    stop("No markers to impute.")
  
  if (n.core > 1) {
    cl <- makeCluster(n.core)
    clusterExport(cl=cl,varlist=NULL)
  }
  
  pred.ix <- matrix(0,nrow=mi,ncol=params$n.mark)
  rownames(pred.ix) <- map.imp$marker
  for (i in 1:mi) {
    ix <- which(map2$chrom==map.imp$chrom[i])
    pred.ix[i,] <- map2$marker[ix[order(abs(map2$pos[ix]-map.imp$pos[i]))[1:params$n.mark]]]
  }
    
  if (n.core > 1) {
    tmp <- parApply(cl,array(map.imp$marker),MARGIN=1,FUN=impute.RF,
                           geno.train=geno1, geno.test=geno2, pred.ix=pred.ix, params=params)
  } else {
    tmp <- apply(array(map.imp$marker),MARGIN=1,FUN=impute.RF,
                        geno.train=geno1, geno.test=geno2, pred.ix=pred.ix, params=params)
  }
  geno.imp <- matrix(sapply(tmp,"[[",1),nrow=mi,byrow = T)
  rownames(geno.imp) <- map.imp$marker
  error <- apply(matrix(sapply(tmp,"[[",2),nrow=mi,byrow = T),2,mean)
  
  if (n.core > 1) {
    stopCluster(cl)
  }
 
  geno.out <- rbind(t(geno2),geno.imp)
  map.out <- rbind(map2,map.imp)
  map.out <- map.out[order(map.out$chrom,map.out$pos),]
  geno.out <- geno.out[map.out$marker,]
  
  write.csv(cbind(map.out,geno.out),file=out.file,row.names=F)
  return(error)
}
