########################
# remove_singularities #
########################

# function to remove the singular variable using a linear model
# before the computation of a mixed model

remove_singularities <- function(d){
  
  fix_form <- paste0('trait~-1 + cross_env+', paste(colnames(d)[5:ncol(d)], collapse = '+'))
  
  m_sg <- lm(as.formula(fix_form), data = d)
  coeff <- coefficients(m_sg)
  if(any(is.na(coeff))){
    d <- d[, -which(colnames(d) %in% names(coeff[is.na(coeff)]))]
  }
  
  return(d)
  
}