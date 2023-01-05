###########
# MM_comp #
###########

# Function to calculate the different mixed models without the QTL effect

# the initial mixed model is:

# y = (Env + cr) + g*e + err

MM_comp <- function(mppData, nEnv, y, cof_mat = NULL, VCOV,
                    maxIter = 100, msMaxIter = 100){
  
  nGeno <- dim(mppData$pheno)[1]
  env <- rep(paste0('E', 1:nEnv), each = nGeno)
  cross <- rep(mppData$cross.ind, nEnv)
  geno <- rep(rownames(mppData$pheno), nEnv)
  cross_env <- paste0(cross, '_', env)
  
  if(!is.null(cof_mat)){
    
    cof_mat <- diag(nEnv) %x% cof_mat
    cof_id <- paste0('cof', 1:ncol(cof_mat))
    colnames(cof_mat) <- cof_id
    
    d <- data.frame(trait = y, env = env, cross_env = cross_env, geno = geno)
    d[, 2:4] <- lapply(d[, 2:4], as.factor)
    d <- data.frame(d, cof_mat)
    
    fix_form <- paste0('trait~-1 + cross_env+', paste(cof_id, collapse = '+'))
    
    # remove the variables that produce singularities
    m_sg <- lm(as.formula(fix_form), data = d)
    coeff <- coefficients(m_sg)
    if(any(is.na(coeff))){
      d <- d[, -which(colnames(d) %in% names(coeff[is.na(coeff)]))]
    }
    cof_id <- paste0('cof', 1:(ncol(d) - 4))
    colnames(d)[5:ncol(d)] <- cof_id
    
    fix_form <- paste0('trait~-1 + cross_env+', paste(cof_id, collapse = '+'))
    
  } else {
    
    d <- data.frame(trait = y, env = env, cross_env = cross_env, geno = geno)
    d[, 2:4] <- lapply(d[, 2:4], as.factor)
    fix_form <- 'trait~cross_env'
    
  }
  
  m <- lme_comp(fix_form = fix_form, VCOV = VCOV, data = d,
                maxIter = maxIter, msMaxIter = msMaxIter)
  
  return(list(model = m, data = d))
  
}