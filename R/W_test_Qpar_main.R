####################
# W_test_Qpar_main #
####################

# Wald test for parental main QTL effect.

W_test_Qpar_main <- function(m, nQTL){
  
  Beta <- m$coefficients$fixed
  Q_ind <- grepl(pattern = 'QTL', x = names(Beta))
  B_QTL <- Beta[Q_ind]
  V_QTL_inv <- qr.solve(m$varFix[Q_ind, Q_ind])
  
  Qeff_nm <- names(B_QTL)
  QTL_id <- paste0('QTL', 1:nQTL)
  
  Qeff_list <- vector(mode = 'list', length = nQTL)
  
  for(i in 1:nQTL){
    
    Q_ind_m_i <- grepl(pattern = QTL_id[i], x = Qeff_nm)
    B_QTL_i <- B_QTL[Q_ind_m_i]
    V_QTL_inv_i <- V_QTL_inv[Q_ind_m_i, Q_ind_m_i]
    par_id <- unlist(lapply(strsplit(x = names(B_QTL_i), split = '_'), `[[`, 3))
    
    # decomposition of individual QTL alleles
    W_Qa <- rep(NA, length(B_QTL_i))
    for(q in 1:length(W_Qa)) W_Qa[q] <- (B_QTL_i[q]^2) * diag(V_QTL_inv_i)[q]
    W_Qa <- pchisq(W_Qa, 1, lower.tail = FALSE)
    names(W_Qa) <- names(B_QTL_i)
    
    Q_tab <- data.frame(par = par_id, Effect = B_QTL_i, p_val = W_Qa,
                        log10P = -log10(W_Qa))
    
    Qeff_list[[i]] <- Q_tab
    
  }
  
  return(Qeff_list)
  
}