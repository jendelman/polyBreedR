#' Impute from low to high density markers
#' 
#' Impute from low to high density markers
#' 
#' Two methods are implemented: Random Forest (RF) and PolyOrigin (PO)
#' 
#' With RF, any missing data are imputed separately for each file, using the population mean (DS) or mode (GT) for each marker. For RF, argument \code{params} recognizes the following parameters: "format","n.tree","n.mark". Format can have values "GT" or "DS" (default = "GT"), n.tree is the number of trees (default = 100), and n.mark is the number of markers to use as predictors (default = 100), chosen based on minimum distance to the target. Classification trees are used for GT and regression trees for DS. Any markers in low.file that are not present in the high.file or with an identical position are discarded.
#' 
#' For method PO, the high density file is the phased parents file from PolyOrigin, which has columns marker, chromosome, position (cM), followed by the phased parental genotypes. Argument \code{params} recognized the following parameters: "ped.file" (default is "ped.csv"), "map.file" (default is NULL), and "format", which can take values "GT" (default) or "AD". The map positions in the high density input file from PO are in cM. To export the imputed data with bp positions, provide name of "map.file" with this information. You must enable command line calls to Julia for the PO method.
#' 
#' @param high.file name of high density file
#' @param low.file name of low density file
#' @param out.file name of CSV output file for imputed data
#' @param method "RF" or "PO" (see Details)
#' @param params list of parameters (see Details)
#' @param n.core multicore processing 
#' 
#' @return for RF method, vector of mean OOB error (vs. number of trees)
#' 
#' @import vcfR
#' @importFrom stats cov sd
#' @importFrom parallel makeCluster stopCluster parApply clusterExport
#' @importFrom randomForest randomForest
#' @importFrom data.table fread fwrite
#' @export

impute_L2H <- function(high.file, low.file, out.file, method, params=list(),
                        n.core=1) {
  
  stopifnot(method %in% c("RF","PO"))
 
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
    prix <- as.character(na.omit(pred.ix[i,]))
    if (params$format=="GT") {
      if (length(unique(y))==1) {
        return(list(pred=rep(y[1],ntest), error=rep(0,params$n.tree)))
      } else {
        ans <- randomForest(y=factor(y), x=geno.train[,prix],
                            xtest=geno.test[,prix],
                            ntree=params$n.tree)
        return(list(pred=as.character(ans$test$predicted), error=ans$err.rate[,1]))
      }
    } else {
      if (sd(y)==0) {
        return(list(pred=rep(y[1],ntest), error=rep(0,params$n.tree)))
      } else {
        ans <- suppressWarnings(randomForest(y=y, x=geno.train[,prix],
                                             xtest=geno.test[,prix],
                                             ntree=params$n.tree))
        return(list(pred=ans$test$predicted, error=ans$mse))
      }
    }
  }
  
  if (method=="RF") {
    np <- names(params)
    if (!("format" %in% np))
      params$format <- "GT"
    if (!("n.tree" %in% np))
      params$n.tree <- 100
    if (!("n.mark" %in% np))
      params$n.mark <- 100
    
    data <- read.vcfR(file=high.file,verbose=FALSE)
    map1 <- data.frame(marker=data@fix[,"ID"], 
                       chrom=data@fix[,"CHROM"], 
                       pos=as.integer(data@fix[,"POS"]))
    map1$marker[which(is.na(map1$marker) | map1$marker==".")] <- paste(map1$chrom,map1$pos,sep="_")
    m1 <- nrow(map1)
    cp1 <- paste(map1$chrom, map1$pos, sep="_")
    
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
    map2$marker[which(is.na(map2$marker) | map2$marker==".")] <- paste(map2$chrom,map2$pos,sep="_")
    
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
    
    pred.ix <- matrix(as.character(NA),nrow=mi,ncol=params$n.mark)
    rownames(pred.ix) <- map.imp$marker
    for (i in 1:mi) {
      ix <- which(map2$chrom==map.imp$chrom[i])
      zz <- min(params$n.mark,length(ix))
      pred.ix[i,1:zz] <- map2$marker[ix[order(abs(map2$pos[ix]-map.imp$pos[i]))[1:zz]]]
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
    
    fwrite(cbind(map.out,geno.out),file=out.file,row.names=F)
    return(error)
  } else {
    #method PO
    np <- names(params)
    if (!("ped.file" %in% np))
      params$ped.file <- "ped.csv"
    if (!("map.file" %in% np))
      params$map.file <- NULL
    if (!("format" %in% np))
      params$format <- "GT"
    
    high <- fread(high.file,check.names=F)
    low <- read.vcfR(low.file,verbose=F)
    meta <- queryMETA(low)
    if (params$format=="GT") {
      GT <- extract.gt(low,element="GT")
      geno <- apply(GT,2,GT2DS)
      rownames(geno) <- rownames(GT)
    } else {
      stopifnot("FORMAT=ID=AD" %in% meta)
      geno <- extract.gt(low,element="AD")
      geno <- gsub(",","|",geno)
    }
    combo <- merge(high,data.frame(marker=rownames(geno),geno,check.names = F),
                   by="marker",all.x=T)
    combo <- combo[order(combo$chromosome,combo$position),]
    if (params$format=="AD")
      combo[is.na(combo)] <- "0|0"
    
    ans <- try(setwd("tmp"),silent=T)
    if (inherits(ans,"try-error")) {
      dir.create("tmp")
    } else {
      setwd("..")
    }
    fwrite(combo,"tmp/PO_geno.csv",na = "NA")
    
    #construct Julia script
    con <- file("tmp/po.jl",open="write")
    writeLines("using PolyOrigin;",con)
    writeLines(sub("X","tmp/PO_geno.csv","genofile=\"X\";"),con)
    writeLines(sub("X",params$ped.file,"pedfile=\"X\";"),con)
    writeLines("polyOrigin(genofile,pedfile,isphysmap=false,refineorder=false,refinemap=false,outstem=\"tmp/imputed\");",con)
    close(con)
    system("julia tmp/po.jl")
    
    imputed <- fread("tmp/imputed_postdoseprob.csv")
    cols <- colnames(geno)
    
    posterior.mean <- apply(imputed[,..cols],2,function(z){
      u <- strsplit(z,split="|",fixed=T)
      v <- sapply(u,function(x){
        x <- as.numeric(x)
        ploidy <- length(x)-1
        round(sum(x*(0:ploidy)),2)
        })
      })
    
    if (!is.null(params$map.file)) {
      map <- fread(params$map.file)[,1:3]
    } else {
      map <- combo[,1:3]
    }
    ix <- match(combo$marker,map$marker)
    
    fwrite(cbind(map[ix,],posterior.mean),out.file)
  }
}