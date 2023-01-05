##################
# cof_mat_reduce #
##################

# Function to reduce the size of the matrix of cofactors. It will only
# keep the cofactors that have a significance below cof_pval_sign

# y: trait value in multiple environments

# cross_mat: within environment cross effect matrix

# cof_mat: cofactor matrix distributed over environments

# cof_pval_sign: minimum significance of the cofactor position
# to be selected.

cof_mat_reduce <- function(y, cross_mat, cof_mat, cof_pval_sign){
  
  # remove singularities first
  X <- cbind(cross_mat, cof_mat)
  m <- lm(y ~ -1 + X)
  X <- X[, !is.na(coefficients(m))]
  
  # Calculate the full model
  m <- lm(y ~ -1 + X)
  
  # detect non-significant cof positions
  cof_pval <- summary(m)$coefficients[, 4]
  cr_id <- grep(pattern = 'cross_env', names(cof_pval))
  X <- X[, -cr_id] # remove the cross term
  cof_pval <- cof_pval[-cr_id]
  
  # select the cofactors
  sel_cof <- cof_pval <= cof_pval_sign
  
  if(any(sel_cof)){
    
    cof_mat <- X[, cof_pval <= cof_pval_sign, drop = FALSE]
    colnames(cof_mat) <- paste0('cof', 1:ncol(cof_mat))
    
  } else {
    
    cof_mat <- NULL
    
    }
  
  return(cof_mat)  
  
}