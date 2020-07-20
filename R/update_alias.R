#' Update names based on alias
#' 
#' Update names based on data frame with alias and preferred name
#'
#' Parameter \code{remove.space} indicates whether blank spaces should be removed before string matching
#'
#' @param x Vector of names to update
#' @param alias Data frame with two columns: first is the preferred name and second is the alias
#' @param remove.space TRUE/FALSE
#' 
#' @return Vector with updated names
#' 
#' @export
update_alias <- function(x,alias,remove.space=TRUE) {
  if (remove.space) {
    x <- gsub(pattern=" ",replacement = "",x=x,fixed=T)
  }
  ix <- which(x %in% alias[,2])
  if (length(ix) > 0) {
    x[ix] <- alias[match(x[ix],alias[,2]),1]
  }
  return(x)
}