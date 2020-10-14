#' Genomic prediction of marker effects from multi-environment trials
#' 
#' Additive and dominance marker effects predicted by BLUP, based on variance components estimated with ASReml-R (must have license)
#' 
#' Best practice for the analysis of datasets comprising multiple experiments follows a two-step approach (Damesa et al. 2017). Step 1 generates a vector of genotype BLUEs (and variance-covariance matrix) for each experiment, taking its design and potentially spatial factors into account. In Step 2, the BLUEs are used as a response variable in a mixed model, and the covariance matrix for the residuals is constrained to equal the direct sum of the covariance matrices of the BLUEs. In the current implementation, the Step 2 model contains three genetic effects: additive, digenic dominant, and independent (to capture higher-order non-additive effects). The digenic dominant term has nonzero mean to allow for heterosis/inbreeding depression (Varona et al. 2018). The model also contains a fixed effect for each experiment and independent random effects for the genotype x environment interaction (an environment may contain multiple experiments). Each element of \code{pheno} corresponds to an experiment and is a list with components named "env" and "y". Env is the name of the environment for that experiment, and "y" is a matrix of dimensions n x (n+1), where n is the number of genotypes in the trial. The first column of the matrix contains the BLUEs, and the remaining n columns are the variance-covariance matrix for the BLUEs. ASReml-R is used to estimate variance components, and BLUPs are computed based on standard formulas (Searle et al. 1992). When \code{D} is not NULL, variances are first estimated for the strictly additive (A) model and used as initial values when fitting the additive + dominance (AD) model. If \code{AIC = TRUE}, then whichever model (A vs. AD) has lower AIC is used for BLUP; otherwise, BLUPs are computed based on the AD model regardless of the AIC values. Individuals in \code{pheno} but not in \code{G} are removed from the analysis; to make predictions for unphenotyped individuals, include them in \code{G} (and \code{D}).
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' @references Varona et al. 2018. Frontiers in Genetics 9: 78. doi:10.3389/fgene.2018.00078
#' @references Searle et al. 1992. Variance Components. doi:10.1002/9780470316856
#' 
#' @param pheno List of matrices containing BLUEs and variance-covariance matrix
#' @param G Returned object from \code{\link{G_mat}} for additive effects
#' @param D Optional, returned object from \code{\link{D_mat}} for dominance effects (default is NULL)
#' @param AIC Boolean variable, whether to select A vs. AD model based on AIC (default is TRUE)
#' @param silent Boolean variable, whether to suppress ASReml-R convergence monitoring (default is FALSE)
#' @param workspace Workspace memory for ASReml-R (default is "128mb")
#' 
#' @return List containing
#' \describe{
#' \item{params}{Matrix with AIC and estimated parameters}
#' \item{indiv}{Matrix of predicted values for the individuals in \code{geno}}
#' \item{markers}{Matrix of predicted effects for the markers in \code{geno}}
#' }
#' 
#' @export
#' @import Matrix
#' @importFrom stats model.matrix

Step2 <- function(pheno,G,D=NULL,AIC=TRUE,silent=FALSE,workspace="128mb") {
  
  requireNamespace("asreml")
  
  id <- rownames(G$mat)
  n.id <- length(id)
  pheno <- lapply(pheno,function(x,id){ix <- which(rownames(x$y) %in% id); x$y <- x$y[ix,c(1,1+ix)]; return(x)},id=id)
  
  Gmat <- G$mat + diag(n.id)*1e-6
  Gfactor <- n.id/(n.id-1)*(mean(diag(Gmat))-mean(Gmat))

  n <- length(pheno)
  blues <- NULL
  Rmat <- vector("list",length=n)
  for (i in 1:n) {
    tmp <- data.frame(expt=i,env=pheno[[i]]$env,id=rownames(pheno[[i]]$y),y=as.numeric(pheno[[i]]$y[,1]),stringsAsFactors = F)
    blues <- rbind(blues,tmp)
    Rmat[[i]] <- pheno[[i]]$y[,-1]
  }
  n <- nrow(blues)
  blues$obs <- factor(paste(blues$id,blues$expt,sep="_"))
  blues$expt <- factor(blues$expt)
  blues$env <- factor(blues$env)
  blues$id <- factor(blues$id,levels=id)
  Rmat <- as.matrix(direct_sum(Rmat))
  attr(Rmat,"INVERSE") <- FALSE
  rownames(Rmat) <- colnames(Rmat) <- as.character(blues$obs)
  
  params <- matrix(0*NA,nrow=3,ncol=9)
  colnames(params) <- c("AIC","Fixed","Vg","Vge","Va","Valpha","Vd","Vbeta","mu_beta")
  rownames(params) <- c("I","I+A","I+A+D")
  
  if (!is.null(D)) {
    Dmat <- D$mat + diag(n.id)*1e-6
    Dfactor <- n.id/(n.id-1)*(mean(diag(Dmat))-mean(Dmat))
    blues$Dsum <- apply(D$coeff,1,sum)[blues$id]
  } else {
    params <- params[1:2,1:6]
  }
  
  print("Baseline model.")
  start.table <- asreml::asreml(data=blues,fixed=y~expt,random=~id+vm(obs,source=Rmat),residual=~idv(units),maxit=30,start.values=T)$vparameters.table
  start.table$Value[2] <- 1
  start.table$Constraint[2] <- "F"
  ans <- asreml::asreml(data=blues,fixed=y~expt,random=~id+vm(obs,source=Rmat),residual=~idv(units),maxit=30,G.param=start.table,na.action=asreml::na.method(y="omit"),trace=!silent,workspace=workspace)
  if (!ans$converge) {
    stop("ASReml-R failed to converge.")
  }
  ans1 <- summary(ans,coef=TRUE)
  params[1,"AIC"] <- round(ans1$aic,1)
  coeff <- ans1$coef.fixed
  params[1,"Fixed"] <- mean(coeff[-match("(Intercept)",rownames(coeff)),1]) + coeff["(Intercept)",1]
  vc1 <- ans1$varcomp
  params[1,c("Vg","Vge")] <- vc1[c(1,3),1]

  #Additive model
  print("Additive model.")
  start.table <- asreml::asreml(data=blues,fixed=y~expt,random=~id+vm(id,source=Gmat,singG="PSD")+vm(obs,source=Rmat),residual=~idv(units),maxit=30,start.values=T)$vparameters.table
  start.table$Value[1] <- vc1$component[1]/2
  start.table$Value[2] <- vc1$component[1]/2
  start.table$Value[3] <- 1
  start.table$Constraint[3] <- "F"
  start.table$Value[5] <- vc1$component[3]
  ans <- asreml::asreml(data=blues,fixed=y~expt,random=~id+vm(id,source=Gmat,singG="PSD")+vm(obs,source=Rmat),residual=~idv(units),maxit=30,G.param=start.table,na.action=asreml::na.method(y="omit"),workspace=workspace,trace=!silent)
  if (!ans$converge) {
    stop("Additive model failed to converge.")
  }
  print("Additive model converged.")
  ans2 <- summary(ans,coef=TRUE)
  params[2,"AIC"] <- round(ans2$aic,1)
  coeff <- ans2$coef.fixed
  params[2,"Fixed"] <- mean(coeff[-match("(Intercept)",rownames(coeff)),1]) + coeff["(Intercept)",1]
  vc2 <- ans2$varcomp
  params[2,c("Va","Vg","Vge","Valpha")] <- vc2[c(2,1,4,2),1]*c(Gfactor,1,1,1/G$scale^2)

  #Dominance model
  if (!is.null(D)) {
    print("Dominance model.")
    start.table <- asreml::asreml(data=blues,fixed=y~expt+Dsum,random=~id+vm(id,source=Dmat,singG="PSD")+vm(id,source=Gmat,singG="PSD")+vm(obs,source=Rmat),residual=~idv(units),maxit=30,start.values=T)$vparameters.table
    start.table$Value[1] <- vc2$component[1]/2
    start.table$Value[2] <- vc2$component[1]/2
    start.table$Value[3] <- vc2$component[2]
    start.table$Value[4] <- 1
    start.table$Constraint[4] <- "F"
    start.table$Value[6] <- vc2$component[4]
    ans <- asreml::asreml(data=blues,fixed=y~expt+Dsum,random=~id+vm(id,source=Dmat,singG="PSD")+vm(id,source=Gmat,singG="PSD")+vm(obs,source=Rmat),residual=~idv(units),maxit=30,G.param=start.table,na.action=asreml::na.method(y="omit"),trace=!silent,workspace=workspace)
    if (!ans$converge) {
      stop("Dominance model failed to converge.")
    }
    print("Dominance model converged.")
    ans3 <- summary(ans,coef=TRUE)
    params[3,"AIC"] <- round(ans3$aic,1)
    coeff <- ans3$coef.fixed
    params[3,"Fixed"] <- mean(coeff[-match(c("Dsum","(Intercept)"),rownames(coeff)),1]) + coeff["(Intercept)",1]  
    vc3 <- ans3$varcomp
    params[3,c("Va","Vd","Vg","Vge","Valpha","Vbeta")] <- 
      vc3[c(3,2,1,5,3,2),1]*c(Gfactor,Dfactor,1,1,1/G$scale^2,1/D$scale^2)
    params[3,"mu_beta"] <- ans3$coef.fixed["Dsum",1]
    
    if (((params[3,"AIC"] < params[2,"AIC"]) & AIC) | !AIC) {
      print("Dominance model used for BLUP.")
      dom.select <- TRUE
    } else {
      dom.select <- FALSE
      print("Additive model used for BLUP.")
      coeff <- ans2$coef.fixed
    }
  }
  
  X <- model.matrix(~expt-1,blues)
  rownames(coeff) <- gsub("_","",rownames(coeff))
  mu_hat <- as.numeric(X %*% coeff[colnames(X),1] + coeff["(Intercept)",1])
  Zg <- model.matrix(~id-1,blues)
  Zge <- sparse.model.matrix(~obs-1,blues)
  if (dom.select) {
    mu_hat <- mu_hat + blues$Dsum*coeff["Dsum",1]
    V <- params[3,"Va"]*tcrossprod(Zg%*%Gmat,Zg) + params[3,"Vd"]*tcrossprod(Zg%*%Dmat,Zg) + params[3,"Vg"]*tcrossprod(Zg) + params[3,"Vge"]*as.matrix(tcrossprod(Zge)) + Rmat
    tmp <- crossprod(Zg,solve(V,blues$y-mu_hat))
    indiv <- data.frame(add=params[3,"Va"]*Gmat%*%tmp,dom=params[3,"Vd"]*Dmat%*%tmp+apply(D$coeff,1,sum)*coeff["Dsum",1])
    markers <- data.frame(add=params[3,"Valpha"]*crossprod(G$coeff,tmp),dom=params[3,"Vbeta"]*crossprod(D$coeff,tmp)+params[3,"mu_beta"])
  } else {
    V <- params[2,"Va"]*tcrossprod(Zg%*%Gmat,Zg) + params[2,"Vg"]*tcrossprod(Zg) + params[2,"Vge"]*as.matrix(tcrossprod(Zge)) + Rmat
    tmp <- crossprod(Zg,solve(V,blues$y-mu_hat))
    indiv <- data.frame(add=params[3,"Va"]*Gmat%*%tmp)
    markers <- data.frame(add=params[3,"Valpha"]*crossprod(G$coeff,tmp))
  }
  return(list(params=params,indiv=indiv,markers=markers))
}