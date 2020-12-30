#' Stage 2 analysis of multi-environment trials (still under development)
#' 
#' Stage 2 analysis of multi-environment trials 
#' 
#' Stage 2 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation (license is required). The variable \code{data} must contain at least three columns: env, id, blue. The first column (env) is the environment identifier, which in plant breeding typically represents a location x year combination. The second column (id) is the genotype identifier, and the third column (blue) is the BLUE from Stage 1 (NAs are not allowed). There are two other reserved column names, which are optional: expt, loc. By default, a fixed effect for each environment is included, but there are situations where BLUEs from multiple experiments (expt) in one environment are included, in which case "expt" overrides "env" to specify the fixed effect portion of the model. When the population of environments includes multiple locations with more than one environment per location, "loc" leads to the inclusion of random effects for genotype x location. For more than 3 locations, a first-order factor-analytic model is used to reduce model complexity. Additional fixed effects can be specified using ASReml-R syntax with the argument \code{fixed} (make sure they have the correct type in \code{data}: numeric vs. factor). To model the uncertainty in the BLUEs from Stage 1 in Stage 2, an additional random effect is included with a constrained variance-covariance matrix named \code{Omega} (following the notation of Damesa et al. 2017). Due to limitations with ASReml-R, this variable must be defined globally instead of passing it to the function. The function \code{\link{Stage2_prep}} can be used to prepare both \code{data} and \code{Omega}. The main effect for genotype can also be partitioned into additive and non-additive effects by defining a global variable named \code{G} for the G matrix. If there are individuals in \code{data} but not in \code{G}, an error is returned. To make predictions for unphenotyped individuals, include them in \code{G}. By default, the workspace and pworkspace limits for ASReml-R are set at 500mb. If you get an error about insufficient memory, try increasing the appropriate value (workspace for variance estimation and pworkspace for BLUP computation).
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' 
#' @param data Data frame with BLUEs from Stage 1 (see Details)
#' @param fixed Additional fixed effects, as a character vector
#' @param workspace Memory limit for ASRreml-R variance estimation
#' @param pworkspace Memory limit for ASRreml-R BLUP computation
#' @param silent TRUE/FALSE, whether to suppress ASReml-R output
#' 
#' @return List containing
#' \describe{
#' \item{aic}{AIC}
#' \item{fixed}{Fixed effect estimates and SE}
#' \item{vc}{Variance component estimates and SE}
#' \item{blup}{BLUPs and reliability (r2) for breeding values if G present or genotypic values if G not present}
#' }
#' 
#' @importFrom stats model.matrix
#' @export

Stage2 <- function(data,fixed=NULL,silent=FALSE,workspace="500mb",pworkspace="500mb") {
  
  stopifnot(requireNamespace("asreml"))
  if (!(exists("Omega")&&inherits(Omega,"matrix"))) {
    stop("Omega covariance matrix is not defined")
  } else {
    stopifnot(nrow(Omega)==nrow(data))
  }
  x <- colnames(data)
  stopifnot(all(c("env","id","blue") %in% x))
  stopifnot(!is.na(data$blue))

  if (exists("G") && inherits(G,"matrix")) {
    GEBV <- TRUE
    print("G matrix detected")
    stopifnot(is.element(as.character(data$id),rownames(G)))
  } else {
    print("No G matrix detected (independent effects)")
    GEBV <- FALSE
  }
  
  if (!is.element("expt",x)) {
    data$expt <- factor(as.character(data$env))
  } else {
    data$expt <- factor(as.character(data$expt))
  }
  data$env <- factor(as.character(data$env))
  data$id <- factor(as.character(data$id))
  Z.id <- model.matrix(~id-1,data)
  colnames(Z.id) <- sub("id","",colnames(Z.id))
  
  #data$env.id <- factor(paste(as.character(data$env),as.character(data$id)))
  X <- model.matrix(~expt-1,data)
  n.fix <- length(fixed)
  if (n.fix > 0) {
    fixed.effects <- paste(c("expt",fixed),collapse="+")
    for (k in 1:n.fix) {
      X <- cbind(X,eval(parse(text=sub("Q",fixed[k],"model.matrix(~Q-1,data)"))))
    }
  } else {
    fixed.effects <- "expt"
  }
  
  random.effects <- "id(id)+vm(units,source=Omega)"
  if (is.element("loc",x)) {
    stop("Not working yet")
    data$loc <- factor(as.character(data$loc))
    Z.gL <- model.matrix(~id:loc-1,data)
    n.loc <- length(levels(data$loc))
    if (n.loc > 3) {
      random.effects <- paste(random.effects,"id(id):fa(loc,1)",sep="+")
    } else {
      random.effects <- paste(random.effects,"id(id):us(loc)",sep="+")
    }
  } 
  if (GEBV) {
    tmp <- paste0(random.effects,"+vm(id,source=G,singG='PSD'),residual=~idv(units),")
  } else {
    tmp <- paste0(random.effects,",residual=~idv(units),")
  }
  
  asreml.options(workspace=workspace,pworkspace=pworkspace,maxit=30,trace=!silent)
  
  print("Estimating variance components...")
  model <- sub(pattern="Q",replacement=fixed.effects,x=paste0("asreml(data=data,fixed=blue~Q,random=~",tmp))
  start.table <- eval(parse(text=paste0(model,"start.values = TRUE)")))$vparameters.table
  k <- grep("Omega",start.table$Component,fixed=T)
  start.table$Value[k] <- 1
  start.table$Constraint[k] <- "F"
  ans <- eval(parse(text=paste0(model,"G.param=start.table)")))
  if (!ans$converge) {
    stop("ASReml-R failed to converge.")
  }
  sans <- summary(ans,coef=TRUE)
  vc <- sans$varcomp
  vc <- vc[vc$bound!="F",c("component","std.error")]
  vc <- rbind(vc,residual=c(mean(diag(Omega)),NA))
  colnames(vc) <- c("estimate","SE")
  rownames(vc) <- replace(rownames(vc),match("units!units",rownames(vc)),"GxE")
  if (GEBV) {
    rownames(vc) <- replace(rownames(vc),grep("vm(id, source = G",rownames(vc),fixed=T),"additive")
    rownames(vc) <- replace(rownames(vc),match("id",rownames(vc)),"non-add")
    k <- match("additive",rownames(vc))
    data$id <- as.character(data$id)
    vc[k,1] <- vc[k,1]*mean(diag(G)[data$id])

    data$id2 <- factor(data$id,levels=rownames(G))
    Z.id2 <- model.matrix(~id2-1,data)
    colnames(Z.id2) <- sub("id2","",colnames(Z.id2))
    Z <- list(Z.id2,Z.id)
    Vu <- list(vc["additive",1]*G,vc["non-add",1]*diag(ncol(Z.id)))
  } else {
    rownames(vc) <- replace(rownames(vc),match("id",rownames(vc)),"G")
    Z <- list(Z.id)
    Vu <- list(vc["G",1]*diag(ncol(Z.id)))
  }
  
  #Prediction using MME
  print("Computing BLUPs...")
  n <- nrow(data)
  ans.MME <- MME(y=data$blue,X=X,Z=Z,Vu=Vu,Rmat=diag(n)*vc["GxE",1]+Omega)
  blup <- ans.MME$blup[[1]]
  blup$BLUP <- blup$BLUP + ans.MME$blue
  fixed <- sans$coef.fixed[,1:2]
  colnames(fixed) <- c("estimate","SE")
  return(list(aic=as.numeric(sans$aic),fixed=fixed,vc=vc,blup=blup))
}