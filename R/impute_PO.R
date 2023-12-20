#' Impute from low to high density markers with PolyOrigin
#' 
#' Impute from low to high density markers by linkage analysis with PolyOrigin
#' 
#' You must have separately installed PolyOrigin and Julia for this function to work.
#'
#'#' The high density file contains phased parental genotypes in PolyOrigin format. The first 3 columns are the genetic map in cM: marker, chrom, position. To output imputed data with physical rather than genetic map positions, including a fourth column named "bp". Subsequent columns are the phased parental genotypes. 
#' 
#' VCF is assumed for the low-density file. The pedigree file must follow PolyOrigin format.
#' 
#' Two output files are generated, corresponding to the posterior maximum and mean genotypes.
#' 
#' @param high.file name of high density file with phased parents
#' @param low.file name of low density VCF file with progeny
#' @param low.format either "GT" (default) or "AD"
#' @param ped.file pedigree file for progeny (must follow PO format)
#' @param out.prefix prefix for output files
#' 
#' @return NULL
#' 
#' @import vcfR
#' @importFrom data.table fwrite
#' @export

impute_PO <- function(high.file, low.file, low.format="GT", ped.file, out.prefix) {

  stopifnot(low.format %in% c("GT","AD"))
  
  high <- fread(high.file,check.names=F)
  low <- read.vcfR(low.file,verbose=F)
  meta <- queryMETA(low)
  if (low.format=="GT") {
    geno <- GT2DS(extract.gt(low,element="GT"))
  } else {
    stopifnot("FORMAT=ID=AD" %in% meta)
    geno <- extract.gt(low,element="AD")
    geno <- gsub(",","|",geno)
  }
  
  if ("bp" %in% colnames(high)) {
    map <- high[,c(1,2,4)]
    high <- high[,-4]
  } else {
    map <- high[,1:3]
  }

  colnames(high) <- replace(colnames(high),1:3,c("marker","chrom","pos"))
  colnames(map) <- c("marker","chrom","pos")
  geno <- geno[rownames(geno) %in% high$marker,] 
    
  combo <- merge(high,data.frame(marker=rownames(geno),geno,check.names = F),
                 by="marker",all.x=T)
  combo <- combo[order(combo$chrom,combo$pos),]
  if (low.format=="AD")
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
  writeLines(sub("X",ped.file,"pedfile=\"X\";"),con)
  writeLines("polyOrigin(genofile,pedfile,isphysmap=false,refineorder=false,refinemap=false,outstem=\"tmp/imputed\");",con)
  close(con)
  system("julia tmp/po.jl")
  
  imputed <- read.csv("tmp/imputed_postdoseprob.csv",check.names=F,as.is=T)
  cols <- colnames(geno)
  
  post.mean <- apply(imputed[,cols],2,function(z){
    u <- strsplit(z,split="|",fixed=T)
    sapply(u,function(x){
      x <- as.numeric(x)
      ploidy <- length(x)-1
      round(sum(x*(0:ploidy)),2)
      })
    })
  fwrite(cbind(map,post.mean),file=sub("Q",out.prefix,"Q_post_mean.csv"))
  
  post.max <- apply(imputed[,cols],2,function(z){
    u <- strsplit(z,split="|",fixed=T)
    sapply(u,function(x){
      which.max(as.numeric(x))-1
    })
  })
  
  fwrite(cbind(map,post.max),file=sub("Q",out.prefix,"Q_post_max.csv"))
}
