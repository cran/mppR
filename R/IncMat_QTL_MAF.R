##################
# IncMat_QTL_MAF #
##################

# function to order the QTL incidence matrix parental and ancestral according
# to MAF per interconnected part. From left to right the MAF increase.

IncMat_QTL_MAF <- function(QTL, Q.eff_i, Q.pos_i, mppData, ref.par = NULL){
  
  if ((Q.eff_i == "par") || (Q.eff_i == "anc")){
    
    if(Q.eff_i == "par"){
      
      # 1. determine the interconnected parts
      
      con.part <- design_connectivity(par_per_cross = mppData$par.per.cross,
                                      plot_des = FALSE)
      
    } else if (Q.eff_i == "anc"){
      
      
      par.clu_i <- mppData$par.clu[Q.pos_i, ]
      par.clu_i <- paste0("A.allele", par.clu_i)
      names(par.clu_i) <- mppData$parents
      
      # change the ref.par into reference allele of the connected part.
      
      if(!is.null(ref.par)){
        
        ref.par <- par.clu_i[ref.par]
        
      }
      
      all.p1 <- par.clu_i[mppData$par.per.cross[, 2]]
      all.p2 <- par.clu_i[mppData$par.per.cross[, 3]]
      
      par.per.cross_i <- cbind(mppData$par.per.cross[, 1], all.p1, all.p2)
      
      con.part <- design_connectivity(par_per_cross = par.per.cross_i,
                                      plot_des = FALSE)
      
      
    }
    
    len.con <- unlist(lapply(X = con.part, FUN = function(x) length(x)))
    con.part <- con.part[names(sort(len.con, decreasing = TRUE))]
    
    # 2. order per MAF within connected parts
    
    Q.mat <- c()
    allele_order <- c()
    con.ind <- c()
    
    for(i in seq_along(con.part)){
      
      # subset results of the connected part
      
      con.part_i <- con.part[[i]]
      QTL_i <- QTL[, con.part_i, drop = FALSE]
      con.ind <- c(con.ind, rep(paste0("c", i), length(con.part_i)))
      
      # order according to number of reference values
      
      all.freq <- apply(X = QTL_i, MARGIN = 2,
                        FUN = function(x) sum(round(x, 3) != 0))
      
      
      if(!is.null(ref.par)){
        
        name.Qmat <- names(sort(all.freq))
        name.Qmat <- name.Qmat[name.Qmat != ref.par]
        name.order <- rev(name.Qmat)
        name.Qmat <- c(name.Qmat, ref.par)
        name.order <- c(ref.par, name.order)
        
        Q.mat <- cbind(Q.mat, QTL_i[, name.Qmat])
        allele_ord_i <- name.order
        
      } else {
        
        name.order <- names(sort(all.freq))
        Q.mat <- cbind(Q.mat, QTL_i[, name.order])
        allele_ord_i <- rev(name.order)
        
      }
      
      
      
      # 3. keep reference and allele order
      
      allele_order <- c(allele_order, allele_ord_i)
      
    }
    
    QTL <- Q.mat
    
    
  } else { allele_order <- con.ind <- NULL}
  
  return(list(QTL = QTL, allele_order = allele_order, con.ind = con.ind))
  
}