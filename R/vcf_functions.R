vcf_prep <- function(vcf.file) {
  
  metadata <- character(1000)
  
  nc <- nchar(vcf.file)
  if (substr(vcf.file,nc-1,nc)==".gz") {
    con <- gzfile(vcf.file,open="r")
  } else {
    con <- file(vcf.file,open="r")
  }
  
  temp <- readLines(con,1)
  h <- 0
  while(substr(temp,1,2)=="##") {
    h <- h+1
    metadata[h] <- temp
    temp <- readLines(con,1)
  }
  header <- temp
  samples <- strsplit(header,split="\t",fixed=TRUE)[[1]][-(1:9)]
  n.sample <- length(samples)
  metadata <- metadata[1:h]
  
  m <- 0
  temp <- readLines(con,1)
  while(!(length(temp)==0 || temp=="")) {
    m <- m + 1
    temp <- readLines(con,1)
  }
  close(con)
  
  ix <- union(grep("##FORMAT",metadata,fixed=T),grep("##INFO",metadata,fixed=T))
  if (length(ix) > 0)
    metadata <- metadata[-ix]
  
  metadata <- append(metadata,
                  c("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">",
                    "##INFO=<ID=DP.AVG,Number=1,Type=Float,Description=\"Average Sample Depth\">",
                    "##INFO=<ID=AF.GT,Number=A,Type=Float,Description=\"Allele Frequency based on GT\">",
                    "##INFO=<ID=AB,Number=1,Type=Float,Description=\"Allelic Bias\">",
                    "##INFO=<ID=SE,Number=1,Type=Integer,Description=\"Sequencing Error (PHRED)\">",
                    "##INFO=<ID=OD,Number=1,Type=Integer,Description=\"OverDispersion (PHRED)\">",
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">",
                    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample Depth\">",
                    "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Posterior Mean Dosage\">",
                     "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
                     header))
                     
  return(list(old.meta=h+1,new.meta=metadata,n.mark=m))
}

vcf_extract <- function(x,field) {
  n <- length(x) - 1
  k <- match(field,strsplit(x[1],split=":",fixed=T)[[1]],nomatch=0)
  if (k==0)
    return(rep(as.numeric(NA),n))
  
  x2 <- strsplit(x[-1],split=":",fixed=T)
  sapply(x2,"[[",k)
}

make_GT <- function(x,ploidy) {
  if (is.na(x)) {
    y <- rep(".",ploidy)
  } else {
    y <- c(rep(0,ploidy-x),rep(1,x))
  }
  return(paste(y,collapse="/"))
}

make_info <- function(x) {
  m <- ncol(x)
  z <- paste0(paste(colnames(x),collapse="=K;"),"=K")
  apply(x,1,function(u,z){
    for (i in 1:m) {
      z <- sub("K",as.character(u[i]),z)
    }
    return(z)
  },z=z)
}
