#' Stage 1 analysis of multi-environment trials
#' 
#' Stage 1 analysis of multi-environment trials
#' 
#' Stage 1 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation (license is required). The variable \code{data} must have one column labeled "id" for the individuals,  one labeled "env" for the environments, plus columns for each of the traits to be analyzed. The data for each environment x trait combination are analyzed independently with a linear mixed model. Argument \code{effects} is a named list of character vectors to specify other effects in the model. Each vector has two elements: the first is "fixed" or "random", and the second is "factor" or "numeric". For example, to include a random block effect, use \code{effects=list(block=c("random","factor"))}. To include stand.count as a numeric covariate, use \code{effects=list(stand.count=c("fixed","numeric"))}. By default, the workspace and pworkspace limits for ASReml-R are set at 500mb. If you get an error about insufficient memory, try increasing the appropriate value (workspace for variance estimation and pworkspace for BLUE computation).
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' 
#' @param data data frame with phenotype data
#' @param traits vector of column names from \code{data}
#' @param effects list of other effects in the model
#' @param workspace memory limit for ASRreml-R variance estimation
#' @param pworkspace memory limit for ASRreml-R BLUE computation
#' @param silent TRUE/FALSE, whether to suppress ASReml-R output
#' 
#' @return List containing
#' \describe{
#' \item{H2}{matrix of broad-sense heritability on a plot basis, for each env x trait combination}
#' \item{blue}{data frame of BLUEs for id x traits}
#' \item{blue.vcov}{list of BLUE + variance-covariance matrices (one matrix per trait)}
#' }
#' 
#' @importFrom stats complete.cases
#' @export

Stage1 <- function(data,traits,effects=NULL,silent=TRUE,workspace="500mb",pworkspace="500mb") {
  
  stopifnot(requireNamespace("asreml"))
  stopifnot(traits %in% colnames(data))
  if (!is.null(effects)) {
    stopifnot(names(effects)%in% colnames(data))
  }
  stopifnot(c("id","env") %in% colnames(data))
  
  data$env <- as.character(data$env)
  data$id <- as.character(data$id)
  envs <- unique(data$env)
  n.env <- length(envs)
  n.trait <- length(traits)
  
  #prepare output variables
  blue.vcov <- vector("list",n.trait)
  names(blue.vcov) <- traits
  blue.all <- NULL
  H2 <- matrix(0,nrow=n.env,ncol=n.trait)
  colnames(H2) <- traits
  rownames(H2) <- envs
  
  #prepare asreml command
  asreml::asreml.options(workspace=workspace,pworkspace=pworkspace,maxit=30,trace=!silent)
  model1 <- "asreml::asreml(data=data1,na.action=asreml::na.method(y='omit',x='omit'),fixed=yyy~F"
  model2 <- "random=~R,residual=~idv(units))"
  
  if (!is.null(effects)) {
    var.type <- sapply(effects,"[",1)
    stopifnot(var.type %in% c("fixed","random"))
    fixed <- names(effects)[var.type=="fixed"]
    random <- names(effects)[var.type=="random"]
    
    var.type <- sapply(effects,"[",2)
    stopifnot(var.type %in% c("factor","numeric"))
    factor.vars <- names(effects)[var.type=="factor"]
    numeric.vars <- names(effects)[var.type=="numeric"]
  } else {
    factor.vars <- numeric.vars <- random <- fixed <- character(0)
  }
  
  n.fixed <- length(fixed)
  n.random <- length(random)
  n.factor <- length(factor.vars)
  n.numeric <- length(numeric.vars)
  
  if (n.fixed > 0) {
    blue.model1 <- sub("F",paste(c(fixed,"id"),collapse="+"),model1,fixed=T)
    blup.model1 <- sub("F",paste(fixed,collapse="+"),model1,fixed=T)
  } else {
    blue.model1 <- sub("F","id",model1,fixed=T)
    blup.model1 <- sub("F","1",model1,fixed=T)
  } 

  if (n.random > 0) {
    blue.model2 <- sub("R",paste(random,collapse="+"),model2,fixed=T)
    blup.model2 <- sub("R",paste(c(random,"idv(id)"),collapse="+"),model2,fixed=T)
  } else {
    blue.model2 <- sub("random=~R,","",model2,fixed=T)
    blup.model2 <- sub("R","idv(id)",model2,fixed=T)
  }
  blue.model <- paste(blue.model1,blue.model2,sep=",")
  blup.model <- paste(blup.model1,blup.model2,sep=",")
  
  for (i in 1:n.trait) {
    cat(sub("X",traits[i],"Trait: X\n"))
    blue <- NULL
    vcov <- vector("list",n.env)
    names(vcov) <- envs
    data$yyy <- data[,traits[i]]
    
    for (j in 1:n.env) {
      cat(sub("X",envs[j],"Env: X\n"))
      ix <- which(data$env==envs[j] & !is.na(data$yyy))
      if (length(ix)==0) {
        cat("Warning: No data present\n")
        H2[j,i] <- NA
      } else {
        data1 <- data[ix,]
        data1$id <- factor(as.character(data1$id))
        if (n.factor > 0) {
          for (q in 1:n.factor) {
            eval(parse(text="data1[,factor.vars[q]] <- factor(as.character(data1[,factor.vars[q]]))"))
          }
        }
        if (n.numeric > 0) {
          for (q in 1:n.numeric) {
            eval(parse(text="data1[,numeric.vars[q]] <- as.numeric(data1[,numeric.vars[q]])"))
          }
        }
        
        blup.ans <- eval(parse(text=blup.model))
        blue.ans <- eval(parse(text=blue.model))
        if (!all(blup.ans$converge,blue.ans$converge)) {
          stop("ASReml-R failed to converge.")
        }
        vc <- summary(blup.ans)$varcomp
        Vg <- vc[match("id!id",rownames(vc)),1]
        Ve <- vc[match("units!units",rownames(vc)),1]
        H2[j,i] <- round(Vg/(Vg+Ve),2)
  
        predans <- asreml::predict.asreml(blue.ans,classify="id",vcov = TRUE)
        vcov[[j]] <- predans$vcov
        tmp <- data.frame(id=as.character(predans$pvals$id),env=envs[j],
                          blue=predans$pvals$predicted.value)
        colnames(tmp) <- c("id","env",traits[i])
        blue <- rbind(blue,tmp)
      }
    }
    
    blue.vcov[[traits[i]]] <- cbind(blue[,traits[i]],direct_sum(vcov[!sapply(vcov,is.null)]))
    rownames(blue.vcov[[traits[i]]]) <- apply(blue[,c("id","env")],1,paste,collapse=":")
    if (is.null(blue.all)) {
      blue.all <- blue
    } else {
      blue.all <- merge(blue,blue.all,all=T)
    }
  }
  blue.all <- blue.all[order(blue.all$env,blue.all$id),]
  return(list(H2=H2,blue=blue.all,blue.vcov=blue.vcov))
}
