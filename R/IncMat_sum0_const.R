#####################
# IncMat_sum0_const #
#####################

# Function to add the constraint to the QTL incidence matrix, cross incidence
# matrix and trait vector for sum to zero constraint.

# mppData : mppData object.

# Q.eff: type of QTL effect.

# Q.list: list of QTL incidence matrices.

# Q.pos: vector of numeric QTL position on the map.

# par.clu: parent clustering matrix.

# cross.mat: cross incidence matrix.

# trait: trait vector.

IncMat_sum0_const <- function(mppData, Q.eff, Q.list, Q.pos, cross.mat, trait){
  
  n_QTL <- length(Q.list)
  
  # count the number of connected parts (number of constraint)
  
  tot_con <- 0
  con_part <- vector(mode = 'list', length = n_QTL)
  len_con_part <- c()
  Qind <- rep(0, n_QTL)
  con.ind <- vector(mode = 'list', length = n_QTL)
  
  for(i in 1:n_QTL){
    
    if(Q.eff == "par"){
      
      con.part <- design_connectivity(par_per_cross = mppData$par.per.cross,
                                      plot_des = FALSE)
      
    } else if (Q.eff == "anc"){
      
      
      par.clu_i <- mppData$par.clu[Q.pos[i], ]
      par.clu_i <- paste0("A.allele", par.clu_i)
      names(par.clu_i) <- mppData$parents
      
      all.p1 <- par.clu_i[mppData$par.per.cross[, 2]]
      all.p2 <- par.clu_i[mppData$par.per.cross[, 3]]
      
      par.per.cross_i <- cbind(mppData$par.per.cross[, 1], all.p1, all.p2)
      
      con.part <- design_connectivity(par_per_cross = par.per.cross_i,
                                      plot_des = FALSE)
      
    }
    
    con_part[[i]] <- con.part 
    
    tot_con <- tot_con + length(con.part)
    
    n_alleles <- unlist(lapply(X = con.part, FUN = length))
    
    len_con_part <- c(len_con_part, n_alleles)
    
    Qind[i] <- sum(n_alleles)
    
    con.ind_k <- c()
    
    for(k in seq_along(con.part)){
      
      con.part_k <- con.part[[k]]
      con.ind_k <- c(con.ind_k, rep(paste0("c", k), length(con.part_k)))
      
    }
    
    con.ind[[i]] <- con.ind_k
    
  }
  
  # for the matrix of constraint
  
  if(tot_con == 1){
    
    const_mat <- list(rep(1, len_con_part))
    
  } else {
    
    const <- rep(paste0('c', 1:tot_con), time = len_con_part)
    const <- factor(x = const, levels = unique(const))
    const_mat <- model.matrix( ~ const - 1, data = const)
    
    Qind <- rep(paste0('Q', 1:n_QTL), time = Qind)
    Qind <- factor(x = Qind, levels = unique(Qind))
    
    const_mat <- split(x = data.frame(const_mat), f = Qind)
    
  }
  
  
  # loop over the QTL incidence matrices to add the constraint to the QTL mat
  
  Q_list_ext <- vector(mode = 'list', length = n_QTL)
  
  for(j in 1:n_QTL){
    
    # order the QTL incidence matrix per connected part
    
    con_part_j <- con_part[[j]]
    con_part_lab <- unlist(con_part_j)
    
    Qmat <- Q.list[[j]]
    Qmat <- Qmat[, con_part_lab]
    
    # add the constraint
    
    const_mat_j <- t(as.matrix(const_mat[[j]]))
    colnames(const_mat_j) <- colnames(Qmat)
    
    Qmat <- rbind(Qmat, const_mat_j)
    
    Q_list_ext[[j]] <- Qmat
    
  }
  
  
  # add the constraint to the cross incidence matrix and the trait.
  
  const_cr_mat <- matrix(0, tot_con, dim(cross.mat)[2])
  colnames(const_cr_mat) <- colnames(cross.mat)
  cross_mat_ext <- rbind(cross.mat, const_cr_mat)
  
  trait_ext <- c(trait, rep(0, tot_con))
  
  return(list(Q.list = Q_list_ext, cross.mat = cross_mat_ext,
              trait = trait_ext, con.ind = con.ind))
  
}
