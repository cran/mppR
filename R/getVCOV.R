###########
# getVCOV #
###########

# get the VCOV stucture from a model calculated with nlme

# parameters:

# mppData: mppData object

# model: results of a mixed model calculated with nlme (lme or gls functions)

# VCOV: character specifiying the VCOV structure

# data: data used to calculatethe mixed model

# nEnv: number of environment

# nGeno: number of genotypes

# NA.rem: remove the rows/columns of the matrix associated with missing (NA)
# values 

# inv: select if you want to get the inverse of the VCOV or not. Default = TRUE

getVCOV <- function(mppData, model, VCOV, data, nEnv, NA.rem = TRUE, inv = TRUE){
  
  N <- dim(data)[1]
  nGeno = dim(mppData$pheno)[1]
  
  if(VCOV == 'CS'){
    
    # diagonal part
    V_e <- (model$sigma)^2
    R <- diag(N) * V_e
    
    # off-diagonal part
    VCOV_est <- getVarCov(model)
    V_g <- VCOV_est[1, 1]
    Z <- model.matrix(~ -1 +  geno, data = data)
    
    VCOV <- (Z %*% t(Z) * V_g) + R
    
  } else if (VCOV == 'CSE') {
    
    # diagonal part
    V_ej <- (c(1.0000000, coef(model$modelStruct$varStruct, unconstrained=F))*model$sigma)^2
    cross_env_comb <- unique(data$cross_env)
    names(V_ej)[1] <- cross_env_comb[!(cross_env_comb %in% names(V_ej))]
    VCOV <- diag(V_ej[data$cross_env])
    
    
  } else if (VCOV == 'CS_CSE'){
    
    # diagonal part
    V_ej <- (c(1.0000000, coef(model$modelStruct$varStruct, unconstrained=F))*model$sigma)^2
    cross_env_comb <- unique(data$cross_env)
    names(V_ej)[1] <- cross_env_comb[!(cross_env_comb %in% names(V_ej))]
    R <- diag(V_ej[data$cross_env])
    
    # off-diagonal part
    VCOV_est <- getVarCov(model)
    V_g <- VCOV_est[1, 1]
    Z <- model.matrix(~ -1 +  geno, data = data)
    
    VCOV <- (Z %*% t(Z) * V_g) + R
    
  } else if (VCOV == 'UN'){
    
    K <- diag(nGeno)
    V_env <- getVarCov(model)
    V_env <- V_env[1:nEnv, 1:nEnv]
    
    VCOV <- V_env %x% K
    VCOV <- VCOV + (diag(N) * ((model$sigma)^2))
    
  }
  
  if(NA.rem){
    
    cc_id <- !is.na(data$trait)
    VCOV <- VCOV[cc_id, cc_id]
    
  }
  
  if(inv){
    
    # VCOV <- solve(VCOV)
    VCOV <- Matrix(VCOV)
    VCOV <- chol2inv(chol(VCOV))
    VCOV <- as.matrix(VCOV)
    
  }
  
  VCOV <- Matrix(VCOV)
  
  return(VCOV)
  
}