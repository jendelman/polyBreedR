#' Stage 2 analysis of multi-environment trials
#' 
#' Stage 2 analysis of multi-environment trials 
#' 
#' Stage 2 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation (license is required). The variable \code{data} must contain at least three columns: id, env, blue. id is the individual identifier, and env represents the environment at which Stage 1 analysis was performed. blue is the BLUE from Stage 1 (NAs are not allowed). Two other column names are reserved: trait and loc. The former triggers a multivariate, multi-trait analysis. The latter triggers the inclusion of a random genotype x location effect (to be completed). To model the uncertainty in the BLUEs from Stage 1 in Stage 2, an additional random effect is included with a variance-covariance matrix named Omega (following notation from Damesa et al. 2017). This variable must be defined globally instead of passing it to the function. The function \code{\link{Stage2_prep}} can be used to prepare both \code{data} and Omega. By default, the model includes independent random effects for genotype (id). Additional genetic effects with specific covariance structure (such as the G matrix for genomic breeding values) can be included using the argument \code{kernels}, which is a vector of variable names (for example, "G") defined in the global environment. (Do not use the name "I" for a kernel; it is reserved for the independent genetic effect.) All individuals in \code{data} must be present in the kernel matrices, but the kernels can contain individuals not in \code{data} to make predictions for unphenotyped individuals using \code{\link{predict_MME}}. All kernel matrices must have the same rownames attribute. For numerical stability when inverting the kernel matrices, a small positive number (1e-5) is added to the diagonal elements. By default, the workspace memory for ASReml-R is set at 500mb. If you get an error about insufficient memory, try increasing it. ASReml-R version 4.1.0.148 or later is required. For kernel matrix K, the variance reported in \code{vars} equal the variance component times the mean of the diagonal elements of ZKZ', which allows for easy computation of the proportion of variance.
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
#' \item{vars}{variances}
#' \item{trait.cov}{genetic variance-covariance matrix for the traits (for multi-trait analysis)}
#' \item{MME}{variable of class \code{\link{MME}} for use with \code{\link{predict_MME}}}
#' }
#' 
#' @importFrom stats model.matrix
#' @importFrom methods new
#' @import Matrix
#' @export

Stage2 <- function(data,kernels=NULL,silent=TRUE,workspace="500mb") {
  
  stopifnot(requireNamespace("asreml"))
  n.obs <- nrow(data)
  if (!(exists("Omega")&&(inherits(Omega,"Matrix")|inherits(Omega,"matrix")))) {
    warning("Omega covariance matrix is not defined")
    skip.omega <- TRUE
    random.effects <- ""
  } else {
    skip.omega <- FALSE
    stopifnot(nrow(Omega)==n.obs)
    random.effects <- "vm(units,source=Omega,singG='PSD')+"
  }
  stopifnot(all(c("id","blue","env") %in% colnames(data)))
  stopifnot(!is.na(data$blue))
  
  if ("trait" %in% colnames(data)) {
    data$Trait <- factor(as.character(data$trait))
    traits <- levels(data$Trait)
    n.trait <- length(traits)
    stopifnot(n.trait>1)
    multi.trait <- TRUE
    fixed.effects <- "env:Trait"
    random.effects <- paste0(random.effects,"id(id):us(Trait)")
    model <- "asreml(data=data,fixed=blue~F,random=~R,residual=~dsum(~id(units)|Trait)"
  } else {
    multi.trait <- FALSE
    fixed.effects <- "env"
    random.effects <- paste0(random.effects,"id(id)")
    model <- "asreml(data=data,fixed=blue~F,random=~R,residual=~idv(units)"
  }
  
  data$env <- factor(as.character(data$env))
  
  nK <- length(kernels) #number of non-identity kernels
  if (nK==0) {
    data$id <- factor(as.character(data$id))
    id <- levels(data$id)
  } else {
    stopifnot(!is.element("I",kernels))
    idK <- vector("list",nK)
    for (i in 1:nK) {
      stopifnot(eval(parse(text=sub("Q",kernels[i],"exists('Q')"))))
      stopifnot(eval(parse(text=gsub("Q",kernels[i],"inherits(Q,'matrix') | inherits(Q,'Matrix')"))))
      idK[[i]] <- eval(parse(text=sub("Q",kernels[i],"rownames(Q)")))
      eval(parse(text=gsub("Q",kernels[i],"colnames(Q) <- rownames(Q)")))
      eval(parse(text=gsub("Q",kernels[i],"Q <- Q + 1e-5*diag(nrow(Q))")))
    }
    stopifnot(sapply(idK,function(z){all(z==idK[[1]])}))
    stopifnot(is.element(levels(data$id),idK[[1]]))
    id <- sort(idK[[1]])
    data$id <- factor(as.character(data$id),levels=id)
    for (i in 1:nK) {
      eval(parse(text=gsub("Q",kernels[i],"Q <- Q[id,id]")))
    }
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
  
  asreml.options(workspace=workspace,maxit=25,trace=!silent)
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
  
  if (multi.trait) {
    vars <- f.id(vc,traits,keyword="!R")
    Rmat <- Diagonal(n=n.obs,x=vars[match(as.character(data$trait),traits)])
    if (!skip.omega) {
      resid <- tapply(diag(Omega),data$Trait,mean)
      vc.out <- rbind(GxE=vars,residual=resid[match(traits,names(resid))])
    } else {
      vc.out <- matrix(vars,nrow=1)
      rownames(vc.out) <- "residual"
    }
    colnames(vc.out) <- traits
  } else {
    vars <- vc[match("units!units",rownames(vc)),1]
    Rmat <- Diagonal(n.obs)*vars
    if (!skip.omega) {
      vc.out <- rbind(GxE=vars,residual=mean(diag(Omega)))
    } else {
      vc.out <- matrix(vars,nrow=1)
      rownames(vc.out) <- "residual"
    }
    colnames(vc.out) <- c("estimate")
  }
  if (!skip.omega) {
    Rmat <- Rmat + Matrix(Omega)
  }
  
  K <- vector("list",nK+1)
  names(K) <- c("I",kernels)
  nZ <- ncol(Z)
  n.id <- length(id)
  
  if (!multi.trait) {
    vcnames <- rownames(vc)
    K[[1]] <- vc["id",1]*Diagonal(n=nZ)
    dimnames(K[[1]]) <- list(id,id)
    vars <- matrix(0,ncol=1,nrow=nK+1)
    rownames(vars) <- paste0("kernel=",c("I",kernels))[(nK+1):1]
    vars[nK+1,1] <- vc["id",1]
    if (nK > 0) {
      for (i in 1:nK) {
        j <- grep(paste("vm(id, source =",kernels[i]),vcnames,fixed=T)
        K[[i+1]] <- eval(parse(text=sub("Q",kernels[i],"Matrix(vc[j,1]*Q)")))
        vars[nK+1-i,1] <- vc[j,1]*eval(parse(text=sub("Q",kernels[i],"mean(diag(Z%*%Q%*%t(Z)))")))
      }
    }
    vc.out <- rbind(vars,vc.out)
  } else {
    trait.cov <- f.cor(vc,traits)
    tmp <- expand.grid(traits,id)
    Knames <- apply(cbind(as.character(tmp$Var2),as.character(tmp$Var1)),1,paste,collapse=":")
    iz <- match(colnames(Z),Knames)
    
    #kernel = I
    if (nK==0) {
      KK <- kronecker(Diagonal(n.id),trait.cov)
      vc.out <- rbind(diag(trait.cov),vc.out)
    } else {
      vars <- f.id(vc=vc,traits=traits,keyword="id:Trait")
      KK <- kronecker(Diagonal(n.id),Diagonal(n=n.trait,x=vars))
      vc.out <- rbind(vars,vc.out)
    }
    dimnames(KK) <- list(Knames,Knames)
    K[[1]] <- KK[iz,iz]
    rownames(vc.out) <- replace(rownames(vc.out),1,"kernel=I")
    
    if (nK > 0) {
      KK <- eval(parse(text=sub("Q",kernels[1],"Matrix(kronecker(Q,trait.cov))")))
      dimnames(KK) <- list(Knames,Knames)
      K[[2]] <- KK[iz,iz]
      tmp <- eval(parse(text=sub("Q","K[[2]]","diag(Z%*%Q%*%t(Z))")))
      vars <- tapply(tmp,data$Trait,mean)
      vc.out <- rbind(vars,vc.out)
      rownames(vc.out) <- replace(rownames(vc.out),1,paste0("kernel=",kernels[1]))
    }
    
    if (nK > 1) {
      for (i in 2:nK) {
        vars <- f.id(vc=vc,traits=traits,keyword=sub("Q",kernels[i],"source = Q"))
        KK <- eval(parse(text=sub("Q",kernels[i],"kronecker(Q,Diagonal(n=n.trait,x=vars))")))
        dimnames(KK) <- list(Knames,Knames)
        K[[i+1]] <- KK[iz,iz]
        tmp <- eval(parse(text=sub("Q","K[[i+1]]","diag(Z%*%Q%*%t(Z))")))
        vars <- tapply(tmp,data$Trait,mean)
        vc.out <- rbind(vars,vc.out)
        rownames(vc.out) <- replace(rownames(vc.out),1,paste0("kernel=",kernels[i]))
      }
    }
  }
  
  if (!multi.trait) {
    out <- new(Class="MME",data=data,kernels=K,Rmat=Rmat)
    return(list(aic=round(as.numeric(sans$aic),1),vars=vc.out,MME=out))
  } else {
    out <- new(Class="MME",data=data[,-match("Trait",colnames(data))],kernels=K,Rmat=Rmat)
    return(list(aic=round(as.numeric(sans$aic),1),vars=vc.out,trait.cov=trait.cov,MME=out))
  }
}
