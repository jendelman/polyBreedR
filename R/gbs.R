#' Genotype calls for GBS 
#' 
#' Genotype calls for genotype-by-sequencing (GBS) data
#' 
#' VCF input file must contain AD field. Posterior mode and mean genotypes are added as GT and DS fields. GQ is also added based on probability of posterior mode. Binomial calculation uses R/updog package (Gerard et al. 2018). Previous INFO is discarded; adds NS, DP.AVG, AF.GT, AB, OD, SE.
#' 
#' @param in.file VCF input file
#' @param out.file VCF output file
#' @param ploidy ploidy
#' @param prior model for prior (see Details)
#' @param bias TRUE/FALSE, whether to estimate allelic bias
#' @param n.core number of cores
#' 
#' @return marker x indiv matrix of read depths
#'
#' @export
#' @importFrom updog flexdog
#' @importFrom parallel makeCluster clusterExport parLapply stopCluster
#' @importFrom stats anova lm chisq.test

gbs <- function(in.file, out.file, ploidy, prior="norm", bias=TRUE, n.core=1) {
  
  if (bias) {
    bias_init <- exp(c(-1, -0.5, 0, 0.5, 1))
  } else {
    bias_init <- 1
  }
  prep <- vcf_prep(in.file)
  
  con.out <- file(out.file,"w")
  for (i in 1:length(prep$new.meta))
    writeLines(con.out,text=prep$new.meta[i])
  
  con.in <- file(in.file,"r")
  tmp <- readLines(con.in,n=prep$old.meta)
  
  if (n.core > 1) {
    cl <- makeCluster(n.core)
    clusterExport(cl=cl,varlist=NULL)
  }
  
  f1 <- function(AD,ploidy,prior) {
    x2 <- strsplit(AD,split=",",fixed=T)
    m <- length(x2)
    ref <- alt <- integer(m)
    ok <- which(sapply(x2,length)==2)
    ref[ok] <- as.integer(sapply(x2[ok],"[[",1))
    alt[ok] <- as.integer(sapply(x2[ok],"[[",2))
    if (any(is.na(ref)))
      ref[is.na(ref)] <- 0
    if (any(is.na(alt)))
      alt[is.na(alt)] <- 0
    
    DP <- ref+alt
    DP.AVG <- paste("DP.AVG",round(mean(DP),1),sep="=")
    
    n <- length(alt)
    tmp <- try(flexdog(refvec=alt,sizevec=alt+ref,
                       ploidy=ploidy,model=prior,bias_init=bias_init,
                       verbose=FALSE,update_bias=bias),silent=TRUE)
    if (!inherits(tmp,"try-error")) {
      AF <- mean(tmp$geno,na.rm=T)/ploidy
      AF.GT <- paste("AF.GT",round(AF,3),sep="=")
      NS <- paste("NS",sum(!is.na(tmp$geno)),sep="=")
      
      if (AF > 0) {
        p <- apply(array(0:ploidy),1,function(z){AF^z*(1-AF)^(ploidy-z)*choose(ploidy,z)})
        #HWE.P <- max(chisq.test(x=table(factor(tmp$geno,levels=0:ploidy)),p=p)$p.value,1e-9)
        AB <- paste("AB",round(tmp$bias,1),sep="=")
        SE <- paste("SE",floor(-10*log10(max(1e-9,tmp$seq))),sep="=")
        OD <- paste("OD",floor(-10*log10(max(1e-9,tmp$od))),sep="=")
        #MIN.DP <- paste("MIN.DP",round(min(tapply(DP,tmp$geno,mean,na.rm=T)),0),sep="=")
        #HWE.P <- paste("HWE.P",round(-10*log10(HWE.P),0),sep="=")
        info <- paste(c(NS,DP.AVG,AF.GT,AB,SE,OD),collapse=";")
      } else {
        info <- paste(c(NS,DP.AVG,AF.GT),collapse=";")
      }
      
      GT <- apply(array(tmp$geno),1,make_GT,ploidy=ploidy)
      DS <- as.character(round(tmp$postmean,1))
      error <- ifelse(tmp$maxpostprob > 1 - 1e-9, 1e-9, 1 - tmp$maxpostprob)
      GQ <- as.character(floor(ifelse(is.na(error),0,-10*log10(error))))
      
    } else {
      GT <- rep(make_GT(as.integer(NA),ploidy=ploidy),n)
      GQ <- DS <- rep(".",n)
      info <- DP.AVG
    }
    z <- apply(cbind(GT,AD,as.character(DP),DS,GQ),1,paste,collapse=":")
    return(paste(c(info,"GT:AD:DP:DS:GQ",z),collapse="\t"))
  }

  block.size <- 100
  nb <- prep$n.mark %/% block.size+1
  i=1
  for (i in 1:nb) {
    tmp <- readLines(con.in,block.size)
    m <- length(tmp)
    cat(sub("X",(i-1)*block.size + m,"Progress: X markers\n"))
    tmp2 <- strsplit(tmp,split="\t",fixed=T)
    x <- lapply(tmp2,function(x){vcf_extract(x[-(1:8)],"AD")})
    if (n.core > 1) {
      ans <- parLapply(cl=cl,X=x,f1,ploidy=ploidy,prior=prior)
    } else {
      ans <- lapply(x,f1,ploidy=ploidy,prior=prior)
    }
    for (j in 1:length(ans)) {
      writeLines(con.out,text=paste(c(tmp2[[j]][1:7],ans[[j]]),collapse="\t"))
    }
  }
  close(con.out)
  close(con.in)
  
  if (n.core > 1)
    stopCluster(cl)
}


