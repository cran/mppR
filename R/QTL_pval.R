############
# QTL_pval #
############

# function to get the QTL decomposed genetic effect for the cross-specific
# parental and ancestral linear models.


QTL_pval <- function(mppData, model, Q.eff, x) {
  
  coeffs <- coef(model)
  index <- which(substr(names(coeffs), 1, 3) == "QTL")
  coeffs <- coeffs[index]
  
  var.comp <- sqrt(diag(vcov(model)))
  index <- which(substr(names(var.comp), 1, 3) == "QTL")
  var.comp <- var.comp[index]
  
  var.comp.full <- rep(NA, length(coeffs))
  var.comp.full[match(names(var.comp), names(coeffs))] <- var.comp
  
  pval <- 2 * pt(q = abs(coeffs/var.comp.full),
                 df = df.residual(model), lower.tail = FALSE)
  pval <- pval * sign(coeffs)
  
  names(pval) <- substr(names(pval), 4, nchar(names(pval)))
  
  if (Q.eff == "cr") {
    
    pval <- pval[paste0("Cr", unique(mppData$cross.ind))]
    
    
  } else if (Q.eff == "par") {
    
    pval <- pval[mppData$parents]
    
  } else if (Q.eff == "anc") {
    
    ref.all <- paste0("A.allele", mppData$par.clu[x, ])
    pval <- pval[ref.all]
    
  }
  
  return(pval)
  
}