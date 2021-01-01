#' S4 class for solving the mixed model equations
#' 
#' @slot y response 
#' @slot X design matrix for the fixed effects
#' @slot Z design matrix for random genetic effects. Colnames must match rownames of the matrices in \code{K}.
#' @slot kernels list of variance-covariance matrices for the genetic effects. Matrices must have the same rownames attribute.
#' @slot Rmat residual variance-covariance matrix

#' 
#' @importFrom methods new
#' @export
MME <- setClass("MME",slots=c(y="numeric",X="Matrix",Z="Matrix",kernels="list",Rmat="Matrix"))