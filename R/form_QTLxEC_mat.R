###################
# form_QTLxEC_mat #
###################

# Function to form the QTLxEC matrix

form_QTLxEC_mat <- function(QTL_list, EC, Q_dist){
  
  nQTL <- length(QTL_list)
  EC_id <- paste0('EC', 1:ncol(EC))
  QTL_EC <- vector(mode = 'list', length = nQTL)
  QTL_EC_nm <- c()
  
  for(i in 1:nQTL){
    
    Qmat_i <- QTL_list[[i]]
    P_sign <- Q_dist[[2]][i, ]
    P_sign <- P_sign[colnames(Qmat_i)]
    
    if(any(P_sign)){
      
      QTL_EC[[i]] <- EC %x% Qmat_i[, P_sign, drop = FALSE]
      
      Q_EC_nm_i <- paste0(paste0('QTL', i), '_', rep(EC_id, each = sum(P_sign)), '_')
      Q_EC_nm_i <- paste0(Q_EC_nm_i, rep(names(P_sign)[P_sign], ncol(EC)))
      QTL_EC_nm <- c(QTL_EC_nm, Q_EC_nm_i)
      
    }
    
  }
  
  QTL_mat_EC <- do.call(cbind, QTL_EC) 
  colnames(QTL_mat_EC) <- QTL_EC_nm
  
  return(QTL_mat_EC)
  
}