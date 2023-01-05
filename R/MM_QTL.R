##########
# MM_QTL #
##########

# function to compute a mixed model for the forward loop. Return model and AIC

MM_QTL  <- function(d, VCOV, maxIter, msMaxIter){
  
  d <- remove_singularities(d)
  Q_id <- colnames(d)[5:ncol(d)]
  fix_form <- paste0('trait~-1 + cross_env+', paste(Q_id, collapse = '+'))
  m <- lme_comp(fix_form = fix_form, VCOV = VCOV, data = d,
                maxIter = maxIter, msMaxIter = msMaxIter)
  
  sum_m <- summary(m)
  AIC_m <- sum_m$AIC
  
  return(list(m = m, AIC = AIC_m))
  
}

