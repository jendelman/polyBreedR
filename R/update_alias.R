#' Update names based on alias
#' 
#' Update names based on data frame with alias and preferred name
#'
#' Parameter \code{remove.space} indicates whether blank spaces should be removed before string matching
#'
#' @param x Vector of names to update
#' @param alias Data frame with two columns: first is the preferred name and second is the alias
#' @param remove.space TRUE/FALSE
#' @param filename update names in file (variable "id")
#' 
#' @return Vector with updated names
#' 
#' @export
#' @importFrom utils write.csv

update_alias <- function(x,alias,remove.space=TRUE,filename=NULL) {
  if (!is.null(filename)) {
    z <- read.csv(filename,check.names=F)
    stopifnot("id" %in% colnames(z))
    x <- z$id
  }
  if (remove.space) {
    x <- gsub(pattern=" ",replacement = "",x=x,fixed=T)
  }
  ix <- which(x %in% alias[,2])
  if (length(ix) > 0) {
    x[ix] <- alias[match(x[ix],alias[,2]),1]
  }
  if (!is.null(filename)) {
    ix <- which(x!=z$id)
    z$id <- x
    write.csv(z,file=filename,row.names=F)
    return(z[ix,])
  } else {
    return(x)
  }
}