#' Merge two genotype matrices and impute missing data (deprecated)
#' 
#' Merge two genotype matrices and impute missing data by BLUP
#' 
#' This function is obsolete. Use \code{\link{impute_L2H}} instead.
#' 
#' Designed to impute from low to high density markers. The BLUP method is equivalent to Eq. 4 of Poland et al. (2012), but this function is not iterative. Additional shrinkage toward the mean is applied if needed to keep the imputed values within the range [0,\code{ploidy}]. Missing data in the input matrices are imputed with the population mean for each marker. If an individual appears in both input matrices, it is renamed with suffixes ".1" and ".2" and treated as two different individuals. Monomorphic markers are removed.
#' 
#' @references Poland et al. (2012) Plant Genome 5:103-113.
#'  
#' @param geno1 Genotype matrix (coded 0...ploidy) with dimensions markers x indiv 
#' @param geno2 Genotype matrix (coded 0...ploidy) with dimensions markers x indiv 
#' @param ploidy Either 2 or 4
#' 
#' @return Imputed genotype matrix (markers x indiv)
#' 
#' @importFrom stats cov sd
#' @export

merge_impute <- function(geno1,geno2,ploidy) {
  
  warning("This function is obsolete and included only for legacy reasons. Use impute_L2H instead.")
  mBLUP <- function(b,W1.L,W1.H_L,V12.L,mu2) {
    W1.H_L[is.na(W1.H_L)] <- 0 
    W1.H <- cbind(W1.L,W1.H_L)
    mu1 <- apply(W1.H,1,mean)
    V11.H <- cov(t(W1.H))
    V11.H_L <- cov(t(W1.H_L))
    Q <- diag(nrow(W1.L))-(1-b)*solve(V11.H,V11.H_L)
    V12.H <- solve(t(Q),b*V12.L)
    return(mu2 + crossprod(V12.H,solve(V11.H,W1.H_L-mu1)))
  }
  
  f1 <- function(x,ploidy) {
    mu <- x[1]
    b1 <- (ploidy-mu)/max(x)
    b2 <- (0-mu)/min(x)
    b <- min(c(1,b1,b2))
    y <- mu+b*x[-1]
    y <- ifelse(y<0,0,y)
    y <- ifelse(y>ploidy,ploidy,y)
    return(y)
  }
  
  id1 <- colnames(geno1)
  n1 <- length(id1)
  id2 <- colnames(geno2)
  n2 <- length(id2)
  id12 <- intersect(id1,id2)
  n12 <- length(id12)
  if (n12 > 0) {
    cat("Warning: Some individuals present in both geno1 and geno2\n")
    id1[match(id12,id1)] <- paste(id12,"1",sep=".")
    colnames(geno1) <- id1
    id2[match(id12,id2)] <- paste(id12,"2",sep=".")
    colnames(geno2) <- id2
  }
  
  mark1 <- rownames(geno1)
  mark2 <- rownames(geno2)
  mark12 <- intersect(mark1,mark2)
  m12 <- length(mark12)
  if (m12==0) {
    stop("There are no markers in common.")
  }
  mark1_2 <- setdiff(mark1,mark2)
  m1_2 <- length(mark1_2)
  mark2_1 <- setdiff(mark2,mark1)
  m2_1 <- length(mark2_1)
  
  geno1a <- rbind(geno1,matrix(NA,nrow=m2_1,ncol=n1))
  rownames(geno1a) <- c(mark1,mark2_1)
  geno2a <- rbind(geno2,matrix(NA,nrow=m1_2,ncol=n2))
  rownames(geno2a) <- c(mark2,mark1_2)
  geno <- cbind(geno1a,geno2a[rownames(geno1a),])
  sdans <- apply(geno,1,sd,na.rm=T)
  monomorphic <- rownames(geno)[which(sdans==0)]
  
  if (length(monomorphic)>0) {
    geno1 <- geno1[!is.element(mark1,monomorphic),]
    geno2 <- geno2[!is.element(mark2,monomorphic),]
  
    mark1 <- rownames(geno1)
    mark2 <- rownames(geno2)
    mark12 <- intersect(mark1,mark2)
    m12 <- length(mark12)
    if (m12==0) {
      stop("There are no markers in common.")
    }
    mark1_2 <- setdiff(mark1,mark2)
    m1_2 <- length(mark1_2)
    mark2_1 <- setdiff(mark2,mark1)
    m2_1 <- length(mark2_1)
  }
  rm(list=c("geno1a","geno2a","geno")) 
  gc() #free up memory
  
  W.L <- scale(t(cbind(geno1[mark12,],geno2[mark12,])),scale=F,center=T)
  W.L[is.na(W.L)] <- 0
  V.L <- cov(t(W.L))
    
  if (m1_2 > 0) {
    ans1_2 <- mBLUP(b = m12/(m12+m1_2),
                    W1.L = W.L[id1,],
                    W1.H_L=scale(t(geno1[mark1_2,]),scale=F,center=T),
                    V12.L=V.L[id1,id2],
                    mu2=apply(W.L[id2,],1,mean))
    marker.means <- apply(geno1[mark1_2,],1,mean,na.rm=T)
    ans1_2 <- apply(rbind(marker.means,ans1_2),2,f1,ploidy=ploidy)
  }
          
  if (m2_1 > 0) {
    ans2_1 <- mBLUP(b = m12/(m12+m2_1),
                    W1.L = W.L[id2,],
                    W1.H_L=scale(t(geno2[mark2_1,]),scale=F,center=T),
                    V12.L=V.L[id2,id1],
                    mu2=apply(W.L[id1,],1,mean))
    marker.means <- apply(geno2[mark2_1,],1,mean,na.rm=T)
    ans2_1 <- apply(rbind(marker.means,ans2_1),2,f1,ploidy=ploidy)
    geno1 <- rbind(geno1,t(ans2_1))
  }
  if (m1_2 > 0) {
    geno2 <- rbind(geno2,t(ans1_2))
  }
  geno <- cbind(geno1,geno2[rownames(geno1),])
  return(impute(geno,method="mean"))
}
