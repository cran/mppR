############
# lme_comp #
############

# function that compute different type of linear mixed models.

lme_comp <- function(fix_form, VCOV, data, maxIter, msMaxIter){
  
  ### CS model
  if (VCOV == 'CS'){
    
    m <- lme(as.formula(fix_form), random = ~ 1 | geno,
             control = list(opt = "optim", maxIter = maxIter, msMaxIter = msMaxIter),
             data = data, na.action = na.omit)
    
    ### CSE model
  } else if (VCOV == 'CSE'){
    
    m <- gls(as.formula(fix_form), weights = varIdent(form = ~ 1 | cross_env),
             control = list(opt = "optim", maxIter = maxIter, msMaxIter = msMaxIter),
                      data = data, na.action = na.omit)
    
    ### CS + CSE model  
  } else if (VCOV == 'CS_CSE'){
    
    m <- lme(as.formula(fix_form), random = ~ 1 | geno,
          weights = varIdent(form = ~ 1 | cross_env),
          control = list(opt = "optim", maxIter = maxIter, msMaxIter = msMaxIter),
          data = data, na.action = na.omit)
    
    ### Unstructured model  
  } else if (VCOV == 'UN'){
    
    m <- lme(as.formula(fix_form), random = list(geno = pdSymm(form = ~ -1 + env)),
            control = list(opt = "optim", maxIter = maxIter, msMaxIter = msMaxIter),
                      data = data, na.action = na.omit)
    
  }
  
  return(m)
  
}