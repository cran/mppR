###################
# form_QEI_QxEC_mat #
###################

# Function to form the QTL design matrix with either QEI term
# or with QTL sensitivity slope given the significance of the QEI term.

form_QEI_QxEC_mat <- function(QTL_list, EC, Q_sign, env_id, ref_env = NULL){
  
  nQTL <- length(QTL_list)
  EC_id <- paste0('EC', 1:ncol(EC))
  QTL_EC <- vector(mode = 'list', length = nQTL)
  QTL_EC_nm <- c()
  nEnv <- nrow(EC)
  
  for(i in 1:nQTL){
    
    Qmat_i <- QTL_list[[i]]
    P_sign <- Q_sign[[1]][i, ]
    P_sign <- P_sign[colnames(Qmat_i)]
    par_nm <- colnames(Qmat_i)
    n_par <- length(par_nm)
    
    # loop over the parents
    QTL_EC_i <- c()
    QTL_nm_i <- c()
    
    for(p in 1:n_par){
      
      if(P_sign[p]){ # QxEC effect
        
        QTL_EC_ip <- EC %x% Qmat_i[, p, drop = FALSE]
        QTL_nm_ip <- paste0(paste0('QTL', i), '_', EC_id, '_')
        QTL_nm_ip <- paste0(QTL_nm_ip, rep(par_nm[p]))
        
      } else { # QEI effect
        
        QTL_EC_ip <- diag(nEnv) %x% Qmat_i[, p, drop = FALSE]
        QTL_nm_ip <- paste0('QTL', i, '_', env_id)
        QTL_nm_ip <- paste0(QTL_nm_ip, '_', rep(par_nm[p], nEnv))
        
        if(!is.null(ref_env)){
          env_ref_pos <- which(env_id == ref_env)
          QTL_EC_ip <- QTL_EC_ip[, -env_ref_pos]
          QTL_nm_ip <- QTL_nm_ip[-env_ref_pos]
        }
        
      }
      
      QTL_EC_i <- cbind(QTL_EC_i, QTL_EC_ip)
      QTL_nm_i <- c(QTL_nm_i, QTL_nm_ip)
      
    }
    
    QTL_EC[[i]] <- QTL_EC_i
    QTL_EC_nm <- c(QTL_EC_nm, QTL_nm_i)
    
  }
  
  QTL_mat <- do.call(cbind, QTL_EC) 
  colnames(QTL_mat) <- QTL_EC_nm
  
  return(QTL_mat)
  
}
