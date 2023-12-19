#' Impute from low to high density markers
#' 
#' Impute from low to high density markers
#' 
#' Two methods are available: Random Forest (RF) and PolyOrigin (PO)
#' 
#' With RF, any missing data are imputed separately for each file, using the population mean (DS) or mode (GT) for each marker. For RF, argument \code{params} recognizes the following parameters: "format","n.tree","n.mark","exclude". \code{format} can have values "GT" or "DS" (default = "GT"). \code{n.tree} is the number of trees (default = 100). \code{n.mark} is the number of markers to use as predictors (default = 100), chosen based on minimum distance to the target. \code{exclude} is a vector of sample names to exclude from the high density file, which is useful for cross-validation. Classification trees are used for GT and regression trees for DS. Both VCF and CSV files are recognized for the input files with method RF, based on the file extension. For CSV files, the first three columns should be marker, chrom, pos, and parameter \code{format} still controls whether classification or regression trees are used. 
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
#' @return for RF method, matrix of OOB error with dimensions markers x trees
#' 
#' @import vcfR
#' @importFrom stats cov sd na.omit
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
  impute.RF <- function(data, geno, train.id, pred.id, params) {
    y <- data$y
    ntest <- nrow(geno) - length(train.id)
    if (params$format=="GT") {
      if (length(unique(y))==1) {
        return(list(pred=rep(y[1],ntest), error=rep(0,params$n.tree)))
      } else {
        ans <- randomForest(y=factor(y), x=geno[train.id,data$pred,drop=FALSE],
                            xtest=geno[pred.id,data$pred,drop=FALSE],
                            ntree=params$n.tree)
        return(list(pred=ans$test$predicted, error=ans$err.rate[,1]))
      }
    } else {
      if (sd(y)==0) {
        return(list(pred=rep(y[1],ntest), error=rep(0,params$n.tree)))
      } else {
        ans <- suppressWarnings(randomForest(y=y, x=geno[train.id,data$pred,drop=FALSE],
                                             xtest=geno[pred.id,data$pred,drop=FALSE],
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
      data[[i]] <- list(y=geno1[,i],
                   pred=map2$marker[ix[order(abs(map2$pos[ix]-map1$pos[i]))[1:zz]]])
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
