#' Trim pedigree
#' 
#' Trim pedigree
#' 
#' Removes founders with a single child, provided the founder is not in the focal population. 
#' 
#' @param ped three column pedigree: id, parent1, parent2
#' @param focal vector of id's for the focal population
#' 
#' @keywords internal
#' @return trimmed pedigree
#' 
trim_ped <- function(ped,focal) {
  
  children.id <- function(id,ped) {
    children <- unique(ped$id[ped$parent1==id | ped$parent2==id])
    return(setdiff(children,0))
  }
  
  founders <- setdiff(ped$id[which(ped$parent1==0 & ped$parent2==0)],focal)
  nf <- length(founders)
  n.child <- numeric(nf)
  for (i in 1:nf) {
    n.child[i] <- length(children.id(id=founders[i],ped=ped))
  } 
  remove.id <- founders[which(n.child==1)]
  if (length(remove.id)==0) {
    return(ped)
  }
  ped2 <- ped[-match(remove.id,ped$id),]
  ped2$parent1[ped2$parent1 %in% remove.id] <- 0
  ped2$parent2[ped2$parent2 %in% remove.id] <- 0
  return(trim_ped(ped2,focal))
}
