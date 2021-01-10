#' Stage 2 analysis of multi-environment trials (still under development)
#' 
#' Stage 2 analysis of multi-environment trials 
#' 
#' Stage 2 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation (license is required). The variable \code{data} must contain at least three columns: env, id, blue. The first column (env) is the environment identifier, which in plant breeding typically represents a location x year combination. The second column (id) is the genotype identifier, and the third column (blue) is the BLUE from Stage 1 (NAs are not allowed). There are two other reserved column names, which are optional: expt, loc. By default, a fixed effect for each environment is included, but there are situations where BLUEs from multiple experiments (expt) in one environment are included, in which case "expt" overrides "env" to specify the fixed effect portion of the model. When the population of environments includes multiple locations with more than one environment per location, "loc" leads to the inclusion of random effects for genotype x location. For more than 3 locations, a first-order factor-analytic model is used to reduce model complexity. Additional fixed effects can be specified using ASReml-R syntax with the argument \code{fixed} (make sure they have the correct type in \code{data}: numeric vs. factor). To model the uncertainty in the BLUEs from Stage 1 in Stage 2, an additional random effect is included with a variance-covariance matrix that must be named "Omega" (notation in Damesa et al. 2017). Due to limitations with ASReml-R, this variable must be defined globally instead of passing it to the function. The function \code{\link{Stage2_prep}} can be used to prepare both \code{data} and Omega. By default, the model includes independent random effects for genotype (id). Additional genetic effects with specific covariance structure (such as the G matrix for genomic breeding values) can be included using the argument \code{kernels}, which is a vector of variable names (for example, "G"). (Do not use the name "I" for a kernel; it is reserved for the independent genetic effect.) All individuals in \code{data} must be present in the kernel matrices, but the kernels can contain individuals not in \code{data} to make predictions for unphenotyped individuals using \code{\link{predict_MME}}. All kernel matrices must have the same rownames attribute. By default, the workspace memory for ASReml-R is set at 500mb. If you get an error about insufficient memory, try increasing it.
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' 
#' @param data Data frame with BLUEs from Stage 1 (see Details)
#' @param kernels Character vector with the names of variance-covariance matrices for genetic effects (see Details)
#' @param fixed Additional fixed effects, as a character vector
#' @param silent TRUE/FALSE, whether to suppress ASReml-R output
#' @param workspace Memory limit for ASRreml-R variance estimation
#' 
#' @return List containing
#' \describe{
#' \item{aic}{AIC}
#' \item{fixed}{Fixed effect estimates and SE}
#' \item{vc}{Variance component estimates and SE}
#' \item{MME}{Variable of class \code{\link{MME}}}
#' }
#' 
#' @importFrom stats model.matrix
#' @importFrom methods new
#' @import Matrix
#' @export

Stage2 <- function(data,fixed=NULL,kernels=NULL,silent=TRUE,workspace="500mb") {
  
  stopifnot(requireNamespace("asreml"))
  if (!(exists("Omega")&&inherits(Omega,"matrix"))) {
    stop("Omega covariance matrix is not defined")
  } else {
    stopifnot(nrow(Omega)==nrow(data))
  }
  x <- colnames(data)
  stopifnot(all(c("env","id","blue") %in% x))
  stopifnot(!is.na(data$blue))
  data$id <- as.character(data$id)
  
  nK <- length(kernels) #number of non-identity kernels
  if (nK > 0) {
    stopifnot(!is.element("I",kernels))
    idK <- vector("list",nK)
    for (i in 1:nK) {
      stopifnot(eval(parse(text=sub("Q",kernels[i],"exists('Q')"))))
      stopifnot(eval(parse(text=sub("Q",kernels[i],"inherits(Q,'matrix')"))))
      idK[[i]] <- eval(parse(text=sub("Q",kernels[i],"rownames(Q)")))
      eval(parse(text=gsub("Q",kernels[i],"colnames(Q) <- rownames(Q)")))
    }
    stopifnot(sapply(idK,function(z){all(z==idK[[1]])}))
    stopifnot(is.element(data$id,idK[[1]]))
    data$id <- factor(data$id,levels=idK[[1]]) #expand iid effect to include all id in kernel matrices
  } else {
    data$id <- factor(data$id)
  }
  id <- levels(data$id)
  Z.id <- Matrix(model.matrix(~id-1,data))
  colnames(Z.id) <- sub("id","",colnames(Z.id),fixed=T)
  
  if (!is.element("expt",x)) {
    data$expt <- factor(as.character(data$env))
  } else {
    data$expt <- factor(as.character(data$expt))
  }
  data$env <- factor(as.character(data$env))
  
  X <- Matrix(model.matrix(~expt-1,data))
  n.fix <- length(fixed)
  if (n.fix > 0) {
    fixed.effects <- paste(c("expt",fixed),collapse="+")
    for (k in 1:n.fix) {
      X <- cbind(X,eval(parse(text=sub("Q",fixed[k],"Matrix(model.matrix(~Q-1,data))"))))
    }
  } else {
    fixed.effects <- "expt"
  }
  
  random.effects <- "id(id)+vm(units,source=Omega,singG='PSD')"
  if (is.element("loc",x)) {
    stop("Not working yet for multiple locations")
    data$loc <- factor(as.character(data$loc))
    Z.gL <- Matrix(model.matrix(~id:loc-1,data))
    n.loc <- length(levels(data$loc))
    if (n.loc > 3) {
      random.effects <- paste(random.effects,"id(id):fa(loc,1)",sep="+")
    } else {
      random.effects <- paste(random.effects,"id(id):us(loc)",sep="+")
    }
  } 
  
  asreml.options(workspace=workspace,maxit=30,trace=!silent)
  cat("Base model (iid genotype effects). Estimating variance components...\n")
  model <- sub(pattern="F",replacement=fixed.effects,"asreml(data=data,fixed=blue~F,random=~R,residual=~idv(units),")
  model1 <- sub(pattern="R",replacement=random.effects,model)
  start.table <- eval(parse(text=paste0(model1,"start.values = TRUE)")))$vparameters.table
  k <- grep("Omega",start.table$Component,fixed=T)
  start.table$Value[k] <- 1
  start.table$Constraint[k] <- "F"
  
  ans <- eval(parse(text=paste0(model1,"G.param=start.table)")))
  if (!ans$converge) {
    stop("ASReml-R failed to converge.")
  }
  sans <- summary(ans,coef=TRUE)
  vc <- sans$varcomp
  
  if (nK > 0) {
    cat("Kernel model. Estimating variance components...\n")
    for (i in 1:nK) {
      random.effects <- paste0(random.effects,sub("Q",kernels[i],"+vm(id,source=Q,singG='PSD')"))
    }
    model2 <- sub(pattern="R",replacement=random.effects,model)
    start.table <- eval(parse(text=paste0(model2,"start.values = TRUE)")))$vparameters.table
    k <- grep("Omega",start.table$Component,fixed=T)
    start.table$Constraint[k] <- "F"
    start.table$Value[match(rownames(vc),start.table$Component)] <- vc$component
    start.table$Value[match("id",start.table$Component)] <- vc["id",1]/(nK+1)
    start.table$Value[grep("singG",start.table$Component,fixed=T)] <- vc["id",1]/(nK+1)
    ans <- eval(parse(text=paste0(model2,"G.param=start.table)")))
    if (!ans$converge) {
      stop("ASReml-R failed to converge.")
    }
    sans <- summary(ans,coef=TRUE)
    vc <- sans$varcomp
  }
  fixed <- sans$coef.fixed[,1:2]
  colnames(fixed) <- c("estimate","SE")
  
  vc <- vc[vc$bound!="F",c("component","std.error")]
  vc <- rbind(vc,residual=c(mean(diag(Omega)),NA))
  colnames(vc) <- c("estimate","SE")
  rownames(vc) <- sub("units!units","GxE",rownames(vc),fixed=T)

  K <- vector("list",nK+1)
  names(K) <- c("I",kernels)
  K[[1]] <- vc["id",1]*Diagonal(length(id))
  colnames(K[[1]]) <- rownames(K[[1]]) <- id
  
  tmp <- rownames(vc)
  if (nK > 0) {
    for (i in 1:nK) {
      j <- grep(paste("vm(id, source =",kernels[i]),tmp,fixed=T)
      tmp <- replace(tmp,j,sub("Q",kernels[i],"kernel=Q"))
      #mean.diag[j] <- eval(parse(text=sub("Q",kernels$name[i],"mean(diag(Q)[data$id])")))
      K[[i+1]] <- eval(parse(text=sub("Q",kernels[i],"Matrix(vc[j,1]*Q)")))
    }
    tmp[match("id",tmp)] <- "kernel=I"
  } else {
    tmp[match("id",tmp)] <- "Genotype"
  }
  rownames(vc) <- tmp
  
  x <- new(Class="MME",y=data$blue,X=X,Z=Z.id,kernels=K,Rmat=Diagonal(nrow(data))*vc["GxE",1]+Matrix(Omega))
  return(list(aic=as.numeric(sans$aic),fixed=fixed,vc=vc,MME=x))
}