#' Stage 2 analysis of multi-environment trials
#' 
#' Stage 2 analysis of multi-environment trials 
#' 
#' Stage 2 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation (license is required). The variable \code{data} has two mandatory column names: id = individual (genotype identifier), and env = environment at which Stage 1 analysis was performed. The argument \code{traits} is a character vector that must match column names in \code{data}. Missing data are allowed in the multi-trait but not the single-trait analysis. For single-trait analysis, an additional random effect can be included to partition the residual and GxE effects. The variance-covariance matrix of this effect must be named Omega (following notation from Damesa et al. 2017) and defined globally in the workspace, rather than passing it to the function (this is due to limitations with ASReml-R). The function \code{\link{Stage2_prep}} can be used to prepare both \code{data} and Omega. By default, the model includes independent random effects for genotype (id). Additional genetic effects with specific covariance structure (such as the G matrix for genomic breeding values) can be included using the argument \code{kernels}, which is a vector of variable names (for example, "G") defined in the global environment. (Do not use the name "I" for a kernel; it is reserved for the independent genetic effect.) All individuals in \code{data} must be present in the kernel matrices, but the kernels can contain individuals not in \code{data} to make predictions for unphenotyped individuals using \code{\link{predict_MME}}. All kernel matrices must have the same rownames attribute. For numerical stability when inverting the kernel matrices, a small positive number (1e-5) is added to the diagonal elements. By default, the workspace memory for ASReml-R is set at 500mb. If you get an error about insufficient memory, try increasing it. ASReml-R version 4.1.0.148 or later is required. For kernel matrix K, the variance reported in \code{vars} equal the variance component times the mean of the diagonal elements of ZKZ', to facilitate proper calculation of the proportion of variance.
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' 
#' @param data data frame of BLUEs from Stage 1 (see Details)
#' @param traits character vector of trait names, matching columns in \code{data}
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
#' @importFrom tidyr pivot_longer
#' @export

Stage2 <- function(data,traits,kernels=NULL,silent=TRUE,workspace="500mb") {
  
  stopifnot(requireNamespace("asreml"))
  stopifnot(length(grep("trait",traits))==0)
  stopifnot(all(c("id","env") %in% colnames(data)))
  
  traits <- sort(traits)
  n.trait <- length(traits)
  if (n.trait > 1) {
    fixed.effects <- "env:trait"
    random.effects <- "id(id):us(trait)"
    model <- sub("blue",paste(traits,collapse=","),
                 "asreml::asreml(data=data,fixed=cbind(blue)~FIXED,random=~RANDOM,residual=~id(units):us(trait)",fixed=T)
  } else {
    stopifnot(!is.na(data[,traits]))
    fixed.effects <- "env"
    random.effects <- "id(id)"
    model <- sub("blue",traits,"asreml(data=data,fixed=blue~FIXED,random=~RANDOM,residual=~idv(units)",fixed=T)
  }
  
  if (!(exists("Omega")&&(inherits(Omega,"Matrix")|inherits(Omega,"matrix"))) | (n.trait > 1)) {
    skip.omega <- TRUE
  } else {
    cat("Omega matrix detected\n")
    skip.omega <- FALSE
    stopifnot(nrow(Omega)==nrow(data))
    random.effects <- paste0(random.effects,"+vm(units,source=Omega,singG='PSD')")
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
    if (n.trait > 1) {
      random.effects <- sub("us(trait)","idh(trait)",random.effects,fixed=T)
      random.effects <- paste0(random.effects,sub("Q",kernels[1],"+vm(id,source=Q,singG='PSD'):us(trait)"))
      if (nK > 1) {
        for (i in 2:nK) {
          random.effects <- paste0(random.effects,sub("Q",kernels[i],"+vm(id,source=Q,singG='PSD'):idh(trait)"))
        }
      }
    } else {
      for (i in 1:nK) {
        random.effects <- paste0(random.effects,sub("Q",kernels[i],"+vm(id,source=Q,singG='PSD')"))
      }
    }
  } 
  
  asreml::asreml.options(workspace=workspace,maxit=30,trace=!silent)
  model <- sub(pattern="FIXED",replacement=fixed.effects,model,fixed=T)
  model <- sub(pattern="RANDOM",replacement=random.effects,model,fixed=T)
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
  
  if (n.trait > 1) {
    trait.cov <- f.cor(vc=vc[grep("units",rownames(vc),fixed=T),],traits)
    Rmat <- kronecker(Diagonal(nrow(data)),trait.cov)
    id.env <- paste(data$id,data$env,sep=":")
    tmp <- expand.grid(traits,id.env)
    Rnames <- apply(cbind(as.character(tmp$Var2),as.character(tmp$Var1)),1,paste,collapse=":")
    
    vars <- diag(trait.cov)
    data2 <- data[,c("id","env",traits)]
    data2 <- pivot_longer(data=data2,cols=match(traits,colnames(data2)),
                          names_to="trait",values_to="blue")
    data2$trait <- factor(data2$trait,levels=traits)
    ix <- which(!is.na(data2$blue))
    data2 <- as.data.frame(data2[ix,])
    tmp <- apply(data2[,c("id","env","trait")],1,paste,collapse=":")
    ix <- match(tmp,Rnames)
    Rmat <- Rmat[ix,ix]
  
    Z <- Matrix(model.matrix(~id:trait-1,data2))
    colnames(Z) <- sub("trait","",colnames(Z),fixed=T)
  } else {
    vars <- vc[match("units!units",rownames(vc)),1]
    data2 <- data[,c("id","env",traits)]
    colnames(data2) <- c("id","env","blue")
    Rmat <- Diagonal(nrow(data2))*vars
    Z <- Matrix(model.matrix(~id-1,data2))
  }
  colnames(Z) <- sub("id","",colnames(Z),fixed=T)
  n.id <- length(id)
  
  if (!skip.omega) {
    vc.out <- rbind(GxE=vars,residual=mean(diag(Omega)))
    Rmat <- Rmat + Matrix(Omega)
  } else {
    vc.out <- matrix(vars,nrow=1)
    rownames(vc.out) <- "residual"
  }
  colnames(vc.out) <- traits

  K <- vector("list",nK+1)
  names(K) <- c("I",kernels)

  if (n.trait==1) {
    vcnames <- rownames(vc)
    K[[1]] <- vc["id",1]*Diagonal(n.id)
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
    tmp <- expand.grid(traits,id)
    Knames <- apply(cbind(as.character(tmp$Var2),as.character(tmp$Var1)),1,paste,collapse=":")
    iz <- match(colnames(Z),Knames)
    
    #kernel = I
    if (nK==0) {
      trait.cov <- f.cor(vc=vc[grep("id:trait",rownames(vc),fixed=T),],traits)
      KK <- kronecker(Diagonal(n.id),trait.cov)
      vc.out <- rbind(diag(trait.cov),vc.out)
    } else {
      vars <- f.id(vc=vc,traits=traits,keyword="id:trait")
      KK <- kronecker(Diagonal(n.id),Diagonal(n=n.trait,x=vars))
      vc.out <- rbind(vars,vc.out)
    }
    dimnames(KK) <- list(Knames,Knames)
    K[[1]] <- KK[iz,iz]
    rownames(vc.out) <- replace(rownames(vc.out),1,"kernel=I")
    
    if (nK > 0) {
      trait.cov <- f.cor(vc=vc[grep("):trait",rownames(vc),fixed=T),],traits)
      KK <- eval(parse(text=sub("Q",kernels[1],"Matrix(kronecker(Q,trait.cov))")))
      dimnames(KK) <- list(Knames,Knames)
      K[[2]] <- KK[iz,iz]
      tmp <- eval(parse(text=sub("Q","K[[2]]","diag(Z%*%Q%*%t(Z))")))
      vars <- tapply(tmp,data2$trait,mean)
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
        vars <- tapply(tmp,data2$trait,mean)
        vc.out <- rbind(vars,vc.out)
        rownames(vc.out) <- replace(rownames(vc.out),1,paste0("kernel=",kernels[i]))
      }
    }
  }
  
  out <- new(Class="MME",data=data2,kernels=K,Rmat=Rmat)
  if (n.trait==1) {
    return(list(aic=round(as.numeric(sans$aic),1),vars=vc.out,MME=out))
  } else {
    return(list(aic=round(as.numeric(sans$aic),1),vars=vc.out,trait.cov=trait.cov,MME=out))
  }
}
