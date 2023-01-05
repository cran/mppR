###################
# W_test_Qpar_GxE #
###################

# Wald test for parental QTLxE effect.

W_test_Qpar_GxE <- function(m, nQTL, nEnv, env_id){
  
  Beta <- m$coefficients$fixed
  Q_ind <- grepl(pattern = 'QTL', x = names(Beta))
  B_QTL <- Beta[Q_ind]
  V_QTL_inv <- qr.solve(m$varFix[Q_ind, Q_ind])
  
  Qeff_nm <- names(B_QTL)
  QTL_id <- paste0('QTL', 1:nQTL)
  
  par_id <- unique(unlist(lapply(strsplit(x = Qeff_nm, split = '_'), `[[`, 3)))
  n_par <- length(par_id)
  
  Qeff_list <- vector(mode = 'list', length = nQTL)
  
  for(i in 1:nQTL){
    
    Q_ind_m_i <- grepl(pattern = QTL_id[i], x = Qeff_nm)
    B_QTL_i <- B_QTL[Q_ind_m_i]
    V_QTL_inv_i <- V_QTL_inv[Q_ind_m_i, Q_ind_m_i]
    
    res_p <- matrix(NA, nrow =  n_par, ncol = (nEnv + 2))
    res_p <- data.frame(par = par_id, res_p)
    
    if(!is.null(env_id)){
      
      colnames(res_p) <- c('par', paste0('Effect_', env_id), 'p_val', 'log10P')
      
    } else {colnames(res_p) <- c('par', paste0('Effect_E', 1:nEnv), 'p_val', 'log10P')}
    
    for(p in 1:n_par){
      nm_mod <- paste0(names(B_QTL_i), '_9v_')
      p_ind <- grepl(pattern = paste0(par_id[p], '_'), x = nm_mod)
      if(any(p_ind)){
        
        B_QTL_ip <- B_QTL_i[p_ind]
        V_QTL_inv_ip <- V_QTL_inv_i[p_ind, p_ind]
        
        # Wald
        W_Q_ip <- t(B_QTL_ip) %*% V_QTL_inv_ip %*% B_QTL_ip
        W_Q_ip <- pchisq(W_Q_ip[1, 1], length(B_QTL_ip), lower.tail = FALSE)
        
        # check all effects are there
        if(!(length(B_QTL_ip) == nEnv)){
          E_id <- paste0('_E', 1:nEnv, '_')
          B_QTL_ip_new <- rep(NA, nEnv)
          for(t in 1:length(E_id)){
            test <- grepl(pattern = E_id[t], names(B_QTL_ip))
            if(any(test)){B_QTL_ip_new[t] <- B_QTL_ip[which(test)]}
          }
          B_QTL_ip <- B_QTL_ip_new
        }
        
        res_p[p, 2:ncol(res_p)] <-  c(B_QTL_ip, W_Q_ip, -log10(W_Q_ip))
        
      }
      
    }
    
    Qeff_list[[i]] <- res_p
    
  }
  
  return(Qeff_list)
  
}