#' Plot G vs. A
#' 
#' Plot marker-based vs. pedigree-based additive relationship coefficients
#' 
#' Useful for finding and correcting pedigree errors. If the G or A coefficient for an individual exceeds the threshold, its name is displayed in the figure. If parentage contains one individual, by default a ggplot2 variable will be returned, but the result can also be written to file. If multiple individuals are present, a filename is required.
#' 
#' @param parentage Data frame of individuals to plot, with 3 columns: id,mother,father
#' @param G Genomic relationship matrix
#' @param A Pedigree relationship matrix
#' @param filename Name of PDF file to save the results (optional for one individual)
#' @param thresh.G Threshold above which names are displayed (default Inf)
#' @param thresh.A Threshold above which names are displayed (default 0.5)
#' 
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices dev.off pdf
#' @export
GvsA <- function(parentage,G,A,filename=NULL,thresh.G=Inf,thresh.A=0.5) {
  
  id <- as.character(parentage[,1])
  n <- length(id)
  if (n > 1 & is.null(filename)) {
    stop("Filename is required.")
  }
  if (!is.null(filename)) {
    filename <- sub(pattern=".pdf",replacement="",filename)
    filename <- sub(pattern=".PDF",replacement="",filename)
    pdf(paste(filename,"pdf",sep="."))
  }
  others <- intersect(colnames(G),colnames(A))
  for (i in 1:n) {
    j <- match(id[i],rownames(A),nomatch=0)
    k <- match(id[i],rownames(G),nomatch=0)
    if (j==0 | k==0) {
      cat(sub("Q",id[i],"Warning: Q is not present in G and A\n"))
    }
    z <- match(id[i],others)
    plotme <- data.frame(id=others[-z],x=A[j,others][-z],y=G[k,others][-z],stringsAsFactors = F)
    plotme$id <- ifelse(plotme$x >= thresh.A | plotme$y >= thresh.G,plotme$id,"")
    ttl <- paste(id[i],paste(as.character(parentage[i,2]),as.character(parentage[i,3]),sep=" / "),sep=" = ")
    p <- ggplot(data=plotme,aes(x=x,y=y)) + geom_point() + geom_text_repel(aes(label=id),show.legend=FALSE,segment.color="red",colour="red") + theme_bw() + ggtitle(ttl) + xlab("A") + ylab("G")
    if (!is.null(filename)){print(p)}
  }
  if (!is.null(filename)) {
    dev.off()
  } else {
    return(p)
  }
}