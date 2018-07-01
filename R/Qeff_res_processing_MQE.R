###########################
# Qeff_res_processing_MQE #
###########################

# Processing of the QTL results for a multi-QTL effect QTL list

# arguments

# model : results of multi-QTL model

# mppData : data object

# Q.list : list of QTLs

# cross.mat : cross intercept incidence matrix

# VCOV : 

# Q.eff

# QTL: list of QTLs

# allele_ref: list of parental (ancestral) alleles used as reference

# con.part: list of connected parts



Qeff_res_processing_MQE <- function(Q.res, Q.eff, Q.pos, con.ind, allele_order,
                                    Q.nb, mppData, VCOV){
  
  
  if (Q.eff == "cr"){
    
    # add additive parent
    
    sign.eff <- sign(Q.res[, 1])
    test.sign <- function(x){ if((x == 0)){NA} else if(x < 0){2} else{3}}
    
    par.add.id <- unlist(lapply(sign.eff, FUN = test.sign))
    Add.parent <- diag(mppData$par.per.cross[1:mppData$n.cr, par.add.id])
    
    Q.res <- data.frame(Q.res, Add.parent, stringsAsFactors = FALSE)
    Q.res[, 1] <- abs(Q.res[, 1])
    
    # add parent genotype (if given)
    
    if(!is.null(mppData$geno.par)){
      
      Par.all <- mppData$geno.par[Q.pos, 5:dim(mppData$geno.par)[2]]
      Par.all <- unlist(Par.all)[Add.parent]
      
      Q.res <- data.frame(Q.res, Par.all, stringsAsFactors = FALSE)
      
    }
    
    
  } else if (Q.eff == "par"){
    
    # reorder the parents to put the reference on the top
    
    par.list.ref <- substr(rownames(Q.res), nchar(as.integer(Q.nb)) + 2 ,
                           nchar(rownames(Q.res)))
    
    index <- match(allele_order, par.list.ref)
    Q.res <- Q.res[index, ]
    
    Q.res <- data.frame(Q.res, Con.part = con.ind, stringsAsFactors = FALSE)
    rownames(Q.res) <- allele_order
    
    # add parent genotype (if given)
    
    if(!is.null(mppData$geno.par)){
      
      Par.all <- mppData$geno.par[Q.pos, 5:dim(mppData$geno.par)[2]]
      
      Par.all <- unlist(Par.all)[allele_order]
      
      Q.res <- data.frame(Q.res, Par.all, stringsAsFactors = FALSE)
      
    }
    
    
  } else if (Q.eff == "anc") {
    
    # order the alleles with the references on top
    
    all.list.ref <- substr(rownames(Q.res), nchar(as.integer(Q.nb)) + 2 ,
                           nchar(rownames(Q.res)))
    
    index <- match(allele_order, all.list.ref)
    Q.res <- Q.res[index, ]
    
    # add connected parts indicator
    
    Q.res <- data.frame(Q.res, Con.part = con.ind, stringsAsFactors = FALSE)
    rownames(Q.res) <- allele_order
    
    
    # project into parents
    
    A.allele <- factor(mppData$par.clu[Q.pos, ])
    A <- model.matrix(~ A.allele - 1)
    A <- A[, rownames(Q.res)]
    
    # separate numeric results from con part information
    
    num.res <- as.matrix(Q.res[, 1:4])
    num.res <- cbind(num.res, 1:dim(num.res)[1])
    
    num.res[is.na(num.res)] <- 9999
    proj.num.res <- A %*% num.res
    proj.num.res[proj.num.res == 9999] <- NA
    rownames(proj.num.res) <- mppData$parents
    
    # add the connected part information
    
    con.part.vec <- Q.res[proj.num.res[, 5], c(5, 6), drop = FALSE]
    
    res <- data.frame(proj.num.res, con.part.vec, stringsAsFactors = FALSE)
    res <- res[order(res[, 5]), ]
    res <- res[, -5]
    
    # add parental score if provided
    
    if(!is.null(mppData$geno.par)){
      
      Par.all <- mppData$geno.par[Q.pos, 5:dim(mppData$geno.par)[2]]
      
      Par.all <- unlist(Par.all)[rownames(res)]
      
      Q.res <- data.frame(res, Par.all, stringsAsFactors = FALSE)
      
    }
    
    
  } else if (Q.eff == "biall"){
    
    
    # project into parents and add genotype score if possible
    
    if(!is.null(mppData$geno.par)){
      
      ref.mat2 <- matrix(rep(c(0, 0, 0, 1), mppData$n.par),
                         nrow = mppData$n.par, byrow = TRUE)
      
      Par.all <- mppData$geno.par[Q.pos, 5:dim(mppData$geno.par)[2]]
      Par.all <- unlist(Par.all)
      
      ref.all <- c(mppData$allele.ref[1, Q.pos, drop = FALSE])
      het.sc <- mppData$allele.ref[c(3, 4), Q.pos]
      
      ind.na <- which(is.na(Par.all))
      ind.ref <- which(Par.all == ref.all)
      ind.het <- which(((Par.all == het.sc[1])|(Par.all == het.sc[2])))
      
      ref.mat2[ind.ref, ] <- matrix(rep(as.matrix(Q.res[, 1:4]),
                                        length(ind.ref)),
                                    nrow = length(ind.ref), byrow = TRUE)
      
      ref.mat2[ind.na, ] <- NA
      ref.mat2[ind.het, 1] <- Q.res[1, 1]/2
      
      # add the sign stars
      
      Sign <- sapply(ref.mat2[, 4], FUN = sign.star)
      ref.mat2 <- data.frame(ref.mat2, Sign, stringsAsFactors = FALSE)
      
      # add parents scores
      
      Q.res <- data.frame(ref.mat2, Par.all, stringsAsFactors = FALSE)
      
      
    } 
    
  }
  
  # add columns names
  
  if(VCOV == "h.err"){
    col.names <- c("Effect", "Std.Err", "t-test", "p-value")
  } else {col.names <- c("Effect", "Std.Err", "W-stat", "p-value")}
  
  colnames(Q.res)[1:4] <- col.names
  
  return(Q.res)
  
}