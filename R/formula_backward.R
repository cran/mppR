####################
# formula_backward #
####################

# function to produce a series of formula for mixed model computation in a
# backward elimination.

# arguments

# Q.names Character vector of QTL position names.

# VCOV Character expression indicating the type of variance covariance structure


formula_backward <- function(Q.names, VCOV){
  
  
  Q.vect <- vector("list", length = length(Q.names))
  
  for(i in seq_along(Q.vect)) Q.vect[[i]] <- c(Q.names[-i], Q.names[i])
  
  if(VCOV == "h.err"){
    
    fbegin <- "trait ~ cross.mat +"
    
  } else if((VCOV == "h.err.as") || (VCOV == "cr.err")){
    
    fbegin <- "trait ~ -1 + cr.mat +"
    
  } else if ((VCOV == "pedigree") || (VCOV == "ped_cr.err")) {
    
    fbegin <- "trait ~ 1 +"
    
  }
  
  
  formula.fct <- function(x, fbegin, VCOV) {
    
    if(VCOV == "h.err"){
      
      paste(fbegin, paste(x, collapse = "+"))
      
    } else {
      
      QTL.el <- vapply(X = x, FUN = function(x) paste0("grp(", x, ")"),
                       character(1))
      paste(fbegin, paste(QTL.el, collapse = "+"))
      
    }
    
    
    
    
  }
  
  lapply(X = Q.vect, formula.fct, fbegin = fbegin, VCOV = VCOV)
  
}