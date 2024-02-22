######################
# determine_QEI_sign #
######################

determine_QEI_sign <- function(Q_sign, par_id, thre_QTL = 2){
  
  res_QxE <- c()
  
  for(q in 1:length(Q_sign)){
    
    Q_sign_q <- Q_sign[[q]]
    sign_QxE <- Q_sign_q$logP_QxE >= thre_QTL
    names(sign_QxE) <- Q_sign_q$par
    sign_QxE <- sign_QxE[par_id]
    
    res_QxE <- rbind(res_QxE, sign_QxE)
    
  }
  
  rownames(res_QxE) <- paste0('Q', 1:nrow(res_QxE))
  
  return(list(QTLxE = res_QxE))
  
}
