#########
# W_QTL #
#########

# Function that calculate the approximate Wald statistics and significnance
# for a multi-allelic QTL position and the same for each QTL alleles.

W_QTL <- function(x, y, Vi, mppData, Q.eff, cross_mat, cof_mat = NULL,
                  nEnv, NA_id, ref_par=NULL){
  
  # form the QTL effect
  QTL <- inc_mat_QTL(x = x, mppData = mppData, Q.eff = Q.eff, order.MAF = TRUE,
                     ref_par = ref_par)
  Q_nm <- colnames(QTL)
  
  # remove singularities in QTL term
  if((Q.eff == 'par') || (Q.eff == 'anc')){
    QTL <- QTL[, -dim(QTL)[2]]
    Q_nm <- Q_nm[-length(Q_nm)]
  }
  
  QTL <- diag(nEnv) %x% QTL
  QTL <- QTL[!NA_id, ]
  QTL_nm <- paste0('Q_', Q_nm)
  QTL_nm <- paste0(rep(QTL_nm, nEnv), rep(paste0('_E', 1:nEnv), each = length(Q_nm)))
  colnames(QTL) <- QTL_nm
  
  if(!is.null(cof_mat)){
    
    X <- cbind(cross_mat, cof_mat, QTL)
    
  } else {
    
    X <- cbind(cross_mat, QTL)
    
  }
  
  # remove singularities
  m_sg <- lm(y ~ -1 + X)
  X <- X[, !is.na(coefficients(m_sg))]
  
  X <- Matrix(X)
  XtX <- t(X) %*% Vi %*% X
  XtX <- as.matrix(XtX)
  V_Beta <- tryCatch(chol2inv(chol(XtX)), error = function(e) NULL)  
  
  if(!is.null(V_Beta)){
    
    V_Beta <- Matrix(V_Beta)
    Beta <- tryCatch(as.matrix(V_Beta %*% t(X) %*% Vi %*% y),
                     error = function(e) NULL)
    
    if(!is.null(Beta)){
      
      Q_ind <- grepl(pattern = 'Q_', x = colnames(X))
      Beta_QTL <- Beta[Q_ind, 1, drop = FALSE]
      Eff_sign <- sign(Beta_QTL[, 1])
      V_QTL <- V_Beta[Q_ind, Q_ind]
      
      V_QTL_inv <- tryCatch(qr.solve(V_QTL), error = function(e) NULL)
      
      if(!is.null(V_QTL_inv)){
        
        # global QTL effect
        W_Q <- t(Beta_QTL) %*% V_QTL_inv %*% Beta_QTL
        pval <- pchisq(W_Q[1, 1], nrow(Beta_QTL), lower.tail = FALSE)
        l_pval <- -log10(pval)
        
        # QTL env effect
        B_nm <- colnames(X)[Q_ind]
        E_id <- paste0('E', 1:nEnv)
        EQ_pval <- rep(NA, nEnv)
        
        for(e in 1:nEnv){
          
          E_ind <- grepl(pattern = E_id[e], x = B_nm)
          Beta_QTL_e <- Beta_QTL[E_ind, , drop = FALSE]
          W_e <- t(Beta_QTL_e) %*% V_QTL_inv[E_ind, E_ind] %*% Beta_QTL_e
          EQ_pval[e] <- pchisq(W_e[1, 1], nrow(Beta_QTL_e), lower.tail = FALSE)
          
        }
        
        names(EQ_pval) <- paste0('Q_pval_E', 1:nEnv)
        
        # decomposition of individual QTL alleles
        W_Qa <- rep(NA, length(Beta_QTL))
        for(i in 1:length(W_Qa)) W_Qa[i] <- (Beta_QTL[i]^2) * diag(V_QTL_inv)[i]
        W_Qa <- pchisq(W_Qa, 1, lower.tail = FALSE)
        W_Qa <- W_Qa * Eff_sign
        names(W_Qa) <- colnames(X)[Q_ind]
        
        # organisation of the results
        
        if(Q.eff == 'cr'){
          
          # not yet
          
        } else if (Q.eff == 'par'){
          
          Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.par)
          par.name <- paste0('Q_', rep(mppData$parents, nEnv), Env_name)
          pvals <- W_Qa[par.name]
          names(pvals) <- par.name
          
        } else if (Q.eff == 'anc'){
          
          ref.all <- paste0("Q_A.allele", mppData$par.clu[x, ])
          Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.par)
          ref.all <- paste0(ref.all, Env_name)
          
          pvals <- W_Qa[ref.all]
          
        }
        
        return(c(l_pval, EQ_pval, pvals))
        
      } else {
        
        if(Q.eff == "cr"){ return(c(0, rep(1, nEnv * (mppData$n.cr + 1))))
          
        } else if (Q.eff == "biall") { return(c(0, rep(1, nEnv)))
          
        } else { return(c(0, rep(1, nEnv * (mppData$n.par + 1)))) }
        
      }
      
      
    } else {
      
      if(Q.eff == "cr"){ return(c(0, rep(1, nEnv * (mppData$n.cr + 1))))
        
      } else if (Q.eff == "biall") { return(c(0, rep(1, nEnv)))
        
      } else { return(c(0, rep(1, nEnv * (mppData$n.par + 1)))) }
      
    }
    
  } else {
    
    if(Q.eff == "cr"){ return(c(0, rep(1, nEnv * (mppData$n.cr + 1))))
      
    } else if (Q.eff == "biall") { return(c(0, rep(1, nEnv)))
      
    } else { return(c(0, rep(1, nEnv * (mppData$n.par + 1)))) }
    
  }
  
}
