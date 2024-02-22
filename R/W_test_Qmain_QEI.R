####################
# W_test_Qmain_QEI #
####################

# Sub-function to determine the significance of main and QTL by environment
# interaction global term.

W_test_Qmain_QEI <- function(m, nQTL, n_par, ref_par){
  
  Beta <- m$coefficients$fixed
  Q_ind <- grepl(pattern = 'QTL', x = names(Beta))
  B_QTL <- Beta[Q_ind]
  V_QTL <- m$varFix[Q_ind, Q_ind]
  
  Qeff_nm <- names(B_QTL)
  QTL_id <- paste0('QTL', 1:nQTL)
  Qeff_list <- vector(mode = 'list', length = nQTL)
  
  for(i in 1:nQTL){
    
    Qeff_list_i <- matrix(NA, nrow = n_par-1, ncol = 2)
    
    # main effect matrix
    Q_ind_m_i <- grepl(pattern = paste0(QTL_id[i], "_main"), x = Qeff_nm)
    B_QTL_m <- B_QTL[Q_ind_m_i]
    V_QTL_m <- V_QTL[Q_ind_m_i, Q_ind_m_i]
    par_id_m <- unlist(lapply(strsplit(x = names(B_QTL_m), split = '_'), `[[`, 3))
    
    # QEI effect matrix
    ind_QEI <- (grepl(pattern = QTL_id[i], x = Qeff_nm)) & !Q_ind_m_i
    B_QTL_QEI <- B_QTL[ind_QEI]
    V_QTL_QEI <- V_QTL[ind_QEI, ind_QEI]
    
    # loop over the parents
    for(q in 1:length(par_id_m)){
      
      WQmq <- (B_QTL_m[q]^2) / V_QTL_m[q, q]
      pval_m <- pchisq(WQmq, 1, lower.tail = FALSE)
      Qeff_list_i[q, 1] <- -log10(pval_m)
      
      # QEI effect sign
      QEI_id_q <- grepl(pattern = par_id_m[q], x = names(B_QTL_QEI))
      B_QTL_QEI_q <- B_QTL_QEI[QEI_id_q]
      V_QTL_QEI_q <- V_QTL_QEI[QEI_id_q, QEI_id_q]
      
      WQEI_q <- t(B_QTL_QEI_q) %*% qr.solve(V_QTL_QEI_q) %*% B_QTL_QEI_q
      pval_QEI <- pchisq(WQEI_q[1, 1], length(B_QTL_QEI_q), lower.tail = FALSE)
      Qeff_list_i[q, 2] <- -log10(pval_QEI)
      
    }
    
    Q_tab <- data.frame(par = c(ref_par[i], par_id_m),
                        logP_main = c(NA, Qeff_list_i[, 1]),
                        logP_QxE = c(NA, Qeff_list_i[, 2]))
    
    Qeff_list[[i]] <- Q_tab
    
  }
  
  return(Qeff_list)
  
}