#' Stage 2 analysis of multi-environment trials
#' 
#' Stage 2 analysis of multi-environment trials 
#' 
#' Stage 2 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation (license is required). The variable \code{data} must contain at least three columns: id, env, blue. id is the individual identifier, and env represents the environment at which Stage 1 analysis was performed. blue is the BLUE from Stage 1 (NAs are not allowed). Two other column names are reserved: trait and loc. The former triggers a multivariate, multi-trait analysis (to be completed). The latter triggers the inclusion of a random genotype x location effect (to be completed). To model the uncertainty in the BLUEs from Stage 1 in Stage 2, an additional random effect is included with a variance-covariance matrix named Omega (following notation from Damesa et al. 2017). This variable must be defined globally instead of passing it to the function. The function \code{\link{Stage2_prep}} can be used to prepare both \code{data} and Omega. By default, the model includes independent random effects for genotype (id). Additional genetic effects with specific covariance structure (such as the G matrix for genomic breeding values) can be included using the argument \code{kernels}, which is a vector of variable names (for example, "G") defined in the global environment. (Do not use the name "I" for a kernel; it is reserved for the independent genetic effect.) All individuals in \code{data} must be present in the kernel matrices, but the kernels can contain individuals not in \code{data} to make predictions for unphenotyped individuals using \code{\link{predict_MME}}. All kernel matrices must have the same rownames attribute. By default, the workspace memory for ASReml-R is set at 500mb. If you get an error about insufficient memory, try increasing it. ASReml-R version 4.1.0.148 or later is required.
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' 
#' @param data data frame of BLUEs from Stage 1 (see Details)
#' @param kernels vector of variable names for variance-covariance matrices of the genetic effects (see Details)
#' @param silent TRUE/FALSE, whether to suppress ASReml-R output
#' @param workspace Memory limit for ASRreml-R variance estimation
#' 
#' @return List containing
#' \describe{
#' \item{aic}{AIC}
#' \item{vc}{Variance component estimates and SE}
#' \item{MME}{Variable of class \code{\link{MME}} for use with \code{\link{predict_MME}}}
#' }
#' 
#' @importFrom stats model.matrix
#' @importFrom methods new
#' @import Matrix
#' @export

Stage2 <- function(data,kernels=NULL,silent=TRUE,workspace="500mb") {
  
  stopifnot(requireNamespace("asreml"))
  if (!(exists("Omega")&&(inherits(Omega,"Matrix")|inherits(Omega,"matrix")))) {
    warning("Omega covariance matrix is not defined")
    skip.omega <- TRUE
    random.effects <- ""
  } else {
    skip.omega <- FALSE
    stopifnot(nrow(Omega)==nrow(data))
    random.effects <- "vm(units,source=Omega,singG='PSD')+"
  }
  stopifnot(all(c("id","blue","env") %in% colnames(data)))
  stopifnot(!is.na(data$blue))
  data$id <- factor(as.character(data$id))
  data$env <- factor(as.character(data$env))
  
  if ("trait" %in% colnames(data)) {
    data$Trait <- factor(as.character(data$trait))
    traits <- levels(data$Trait)
    n.trait <- length(traits)
    stopifnot(n.trait>1)
    multi.trait <- TRUE
    fixed.effects <- "env:Trait"
    random.effects <- paste0(random.effects,"id(id):us(Trait)")
    model <- "asreml(data=data,fixed=blue~F,random=~R,residual=~dsum(~id(units)|Trait)"
    X <- Matrix(model.matrix(~env:Trait-1,data))
  } else {
    multi.trait <- FALSE
    fixed.effects <- "env"
    random.effects <- paste0(random.effects,"id(id)")
    model <- "asreml(data=data,fixed=blue~F,random=~R,residual=~idv(units)"
    X <- Matrix(model.matrix(~env-1,data))
  }
  
  nK <- length(kernels) #number of non-identity kernels
  if (nK > 0) {
    data$id <- as.character(data$id)
    stopifnot(!is.element("I",kernels))
    idK <- vector("list",nK)
    for (i in 1:nK) {
      stopifnot(eval(parse(text=sub("Q",kernels[i],"exists('Q')"))))
      stopifnot(eval(parse(text=gsub("Q",kernels[i],"inherits(Q,'matrix') | inherits(Q,'Matrix')"))))
      idK[[i]] <- eval(parse(text=sub("Q",kernels[i],"rownames(Q)")))
      eval(parse(text=gsub("Q",kernels[i],"colnames(Q) <- rownames(Q)")))
    }
    stopifnot(sapply(idK,function(z){all(z==idK[[1]])}))
    stopifnot(is.element(data$id,idK[[1]]))
    data$id <- factor(data$id,levels=idK[[1]]) #expand iid effect to include all id in kernel matrices
    
    if (multi.trait) {
      random.effects <- sub("us(Trait)","idh(Trait)",random.effects,fixed=T)
      random.effects <- paste0(random.effects,sub("Q",kernels[1],"+vm(id,source=Q,singG='PSD'):us(Trait)"))
      if (nK > 1) {
        for (i in 2:nK) {
          random.effects <- paste0(random.effects,sub("Q",kernels[i],"+vm(id,source=Q,singG='PSD'):idh(Trait)"))
        }
      }
    } else {
      for (i in 1:nK) {
        random.effects <- paste0(random.effects,sub("Q",kernels[i],"+vm(id,source=Q,singG='PSD')"))
      }
    }
  } 
  
  if (multi.trait) {
    Z <- Matrix(model.matrix(~id:Trait-1,data))
    colnames(Z) <- sub("Trait","",colnames(Z),fixed=T)
  } else {
    Z <- Matrix(model.matrix(~id-1,data))
  }
  colnames(Z) <- sub("id","",colnames(Z),fixed=T)
  
  asreml.options(workspace=workspace,maxit=30,trace=!silent)
  model <- sub(pattern="F",replacement=fixed.effects,model)
  model <- sub(pattern="R",replacement=random.effects,model)
  if (!skip.omega) {
    start.table <- eval(parse(text=paste0(model,",start.values = TRUE)")))$vparameters.table
    k <- grep("Omega",start.table$Component,fixed=T)
    start.table$Value[k] <- 1
    start.table$Constraint[k] <- "F"
    ans <- eval(parse(text=paste0(model,",G.param=start.table)")))
  } else {
    ans <- eval(parse(text=paste0(model,")")))
  }
  if (!ans$converge) {
    stop("ASReml-R failed to converge.")
  }
  sans <- summary(ans,coef=TRUE)
  vc <- sans$varcomp
  vc <- vc[vc$bound!="F",c("component","std.error")]
  colnames(vc) <- c("estimate","SE")
  if (multi.trait) {
    if (!skip.omega) {
      rownames(vc) <- sub("!R","_GxE",rownames(vc),fixed=T)
    } else {
      rownames(vc) <- sub("!R","_residual",rownames(vc),fixed=T)
    }
  } else {
    if (!skip.omega) {
      rownames(vc) <- sub("units!units","GxE",rownames(vc),fixed=T)
    } else {
      rownames(vc) <- sub("units!units","residual",rownames(vc),fixed=T)
    }
  }
  if (!skip.omega) {
    vc <- rbind(vc,residual=c(mean(diag(Omega)),NA))
    Rmat <- Diagonal(nrow(data))*vc["GxE",1]+Matrix(Omega)
  } else {
    Rmat <- Diagonal(nrow(data))*vc["residual",1]
  }
  
  K <- vector("list",nK+1)
  names(K) <- c("I",kernels)
  vcnames <- rownames(vc)
  if (!multi.trait) {
    K[[1]] <- vc["id",1]*Diagonal(n=ncol(Z))
    dimnames(K[[1]]) <- list(colnames(Z),colnames(Z))
    if (nK > 0) {
      for (i in 1:nK) {
        j <- grep(paste("vm(id, source =",kernels[i]),vcnames,fixed=T)
        vcnames <- replace(vcnames,j,sub("Q",kernels[i],"kernel=Q"))
        K[[i+1]] <- eval(parse(text=sub("Q",kernels[i],"Matrix(vc[j,1]*Q)")))
      }
    }
    vcnames[match("id",vcnames)] <- "kernel=I"
    rownames(vc) <- vcnames
  } else {
    if (nK > 0) {
      ix <- match(apply(array(traits),1,function(trait){sub("Q",trait,"id:Trait!Trait_Q")}),rownames(vc))
      tmp <- lapply(as.list(traits),grep,x=colnames(Z),fixed=T)
      K[[1]] <- Diagonal(n=ncol(Z),x=rep(vc[ix,1],times=sapply(tmp,length)))
      dimnames(K[[1]]) <- list(colnames(Z),colnames(Z))
      
      # first kernel has us(Trait)
      j <- grep(paste("vm(id, source =",kernels[1]),vcnames,fixed=T)  #
      trait_cov <- matrix(0,n.trait,n.trait)
      dimnames(trait_cov) <- list(traits,traits)
      tmp <- sapply(strsplit(vcnames[j],split="Trait!Trait_",fixed=T),"[",2)
      tmp <- strsplit(tmp,split=":",fixed=T)
      trait_cov[cbind(sapply(tmp,"[",1),sapply(tmp,"[",2))] <- vc[j,1]
      trait_cov[cbind(sapply(tmp,"[",2),sapply(tmp,"[",1))] <- vc[j,1]
      K[[2]] <- eval(parse(text=sub("Q",kernels[1],"Matrix(kronecker(Q,trait_cov))")))
      tmp <- expand.grid(x=idK[[1]],y=traits)
      tmp2 <- apply(tmp,1,paste,collapse=":")
      dimnames(K[[2]]) <- list(tmp2,tmp2)
      
      #more kernels FIX THIS
      if (nK > 1) {
        for (i in 2:nK) {
        }
      }
    }
  }
  out <- new(Class="MME",y=data$blue,X=X,Z=Z,kernels=K,Rmat=Rmat)
  return(list(aic=round(as.numeric(sans$aic),1),vc=vc,MME=out))
}