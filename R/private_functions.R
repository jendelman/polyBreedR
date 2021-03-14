f.id <- function(vc,traits,keyword) {
  vcnames <- rownames(vc)
  ix1 <- grep(keyword,vcnames,fixed=T)
  ix2 <- lapply(as.list(traits),grep,x=vcnames,fixed=T)
  ix <- sapply(ix2,intersect,y=ix1)
  return(vc[ix,1])
}

f.cor <- function(vc,traits) {
  tmp <- expand.grid(traits,traits)
  tmp <- apply(tmp,2,as.character)[,c(2,1)]
  tmp <- tmp[tmp[,1] >= tmp[,2],]
  tmp2 <- apply(tmp,1,paste,collapse=":")
  ix <- sapply(as.list(tmp2),grep,x=rownames(vc),fixed=T)
  n.trait <- length(traits)
  trait.cov <- matrix(0,nrow=n.trait,n.trait)
  dimnames(trait.cov) <- list(traits,traits)
  trait.cov[cbind(tmp[,1],tmp[,2])] <- vc[ix,1]
  trait.cov[upper.tri(trait.cov,diag=F)] <- trait.cov[lower.tri(trait.cov,diag=F)]
  return(trait.cov)
}