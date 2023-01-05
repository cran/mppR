############################
# determine_Qmain_QxE_sign #
############################

# function to determine if the parental effect are significant as main effect
# and as QTLxE effect.

determine_Qmain_QxE_sign <- function(Qeff, par_id, thre_QTL = 2, QxE_only = FALSE){
  
  res_main <- c()
  res_QxE <- c()
  
  for(q in 1:length(Qeff)){
    
    Qeff_q <- Qeff[[q]]
    
    
    sign_main <- Qeff_q$logP_main >= thre_QTL
    sign_QxE <- Qeff_q$logP_QxE >= thre_QTL
    QxE_more_main <- Qeff_q$logP_QxE > Qeff_q$logP_main
    
    sign_a <- sign_main | sign_QxE
    
    if(QxE_only){
      
      sign_a_QxE <- sign_QxE
      
    } else {
      
      sign_a_QxE <- sign_a & QxE_more_main
      
    }
    
    names(sign_a) <- names(sign_a_QxE) <- par_id
    
    res_main <- rbind(res_main, sign_a)
    res_QxE <- rbind(res_QxE, sign_a_QxE)
    
  }
  
  rownames(res_main) <- rownames(res_QxE) <- paste0('Q', 1:nrow(res_main))
  
  return(list(QTL_main = res_main, QTLxE = res_QxE))
  
}