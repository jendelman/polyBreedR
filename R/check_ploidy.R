#' Check ploidy
#' 
#' Fraction of simplex or triplex markers 
#' 
#' For every indiv in the genotype matrix, the fraction of markers per chromosome called as simplex or triplex is calculated, which should be low for diploids. A small amount of missing genotype data can be tolerated.
#' 
#' @param geno Genotype matrix (markers x indiv)
#' @param map Data frame with marker map (Marker, Chrom, Position)
#'
#' @return List containing
#' \describe{
#' \item{mat}{Matrix (indiv x chrom) of results}
#' \item{plot}{ggplot2 barplot}
#' }
#'
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export
#' 
check_ploidy <- function(geno,map) {
  colnames(map) <- c("Marker","Chrom","Position")
  marks <- intersect(map$Marker,rownames(geno))
  map <- map[which(map$Marker %in% marks),]
  geno <- geno[map$Marker,]
  
  x <- split(map$Marker,map$Chrom)
  ans <- lapply(x,function(y){apply(geno[y,],2,function(z){ tab <- table(factor(z,levels=0:4))
                                                            sum(tab[c(2,4)])/sum(tab)
                                                           })  })
  ans <- matrix(unlist(ans),nrow=ncol(geno),ncol=length(x),byrow=F)
  rownames(ans) <- colnames(geno)
  colnames(ans) <- names(x)
  plotme <- pivot_longer(data.frame(id=rownames(ans),ans),cols=1+1:ncol(ans),names_to="Chrom",values_to="y")
  plotme$Chrom <- factor(plotme$Chrom)
  v <- apply(ans,1,sum)
  plotme$id <- factor(plotme$id,levels=rownames(ans)[order(v)],ordered=T)
  p <- ggplot(data=plotme,aes(x=id,y=y,fill=Chrom)) + geom_col() + scale_fill_brewer(palette="Set3") + theme_bw() + theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) + ylab("Proportion of Markers with Dosage 1 or 3")
  return(list(mat=ans,plot=p))
}