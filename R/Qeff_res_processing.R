#######################
# Qeff_res_processing #
#######################

# Processing of the QTL results

# arguments

# model : multi-QTL model

# mppData : data object

# Q.list : list of QTLs

# cross.mat : cross intercept incidence matrix

# VCOV : 

# Q.eff

# QTL: list of QTLs

# allele_ref: list of parental (ancestral) alleles used as reference

# con.part: list of connected parts

Qeff_res_processing <- function(model, mppData, cross.mat, Q.list, QTL,
                                Q.eff, VCOV, allele_order, con.ind){
  
  n.QTL <- length(Q.list)
  
  if(is.character(QTL)){ QTL.list <- QTL} else {QTL.list <- QTL[, 1]}
  
  if(VCOV == "h.err"){
    
    results <- summary(model)$coefficients
    index <- (substr(rownames(results), 1, 1) == "Q")
    results <- subset(x = results, subset = index, drop = FALSE)
    
    # possibility to define here the reference names using names(coef(model))
    # Then make the filling of the NA values.
    
    # Then we will assume that for the mixed model we get automatically
    # all estimated values...
    
  } else {
    
    # index <- substr(names(rev(model$coefficients$fixed)), 1, 1) == "Q"
    # 
    # w.table <- asreml::wald(model)
    # w.stat <- w.table[substr(rownames(w.table), 1, 1) == "Q", c(3, 4)]
    # 
    # results <- cbind(rev(model$coefficients$fixed)[index],
    #                  rev(sqrt(model$vcoeff$fixed))[index], w.stat)
    # results <- as.matrix(results)
    
    if(Q.eff == "biall"){ rownames(results) <- paste0("Q", 1:n.QTL) }
    
  }
  
  # control for singular values and fill missing values
  
  if (Q.eff == "cr"){
    
    ref.mat <- matrix(rep(c(0, 0, 0, 1), (mppData$n.cr * n.QTL)),
                      nrow = (mppData$n.cr * n.QTL), byrow = TRUE)
    ref.names <- paste0(rep(paste0("Q", 1:n.QTL), each = mppData$n.cr),
                        rep(colnames(cross.mat), times = n.QTL))
    
    index <- match(rownames(results), ref.names)
    ref.mat[index, ] <- results
    rownames(ref.mat) <- ref.names
    
    # add sign stars
    
    Sign <- sapply(ref.mat[, 4], FUN = sign.star)
    
    results <- data.frame(ref.mat, Sign, stringsAsFactors = FALSE)
    
    Qeff.mat <- vector(mode = "list", n.QTL)
    
    Q.id <- paste0("Q", 1:n.QTL)
    Q.ind <- rep(Q.id, each = mppData$n.cr)
    
    for(i in 1:n.QTL){
      
      # subset QTL
      
      Qi <- results[Q.ind == Q.id[i], ]
      
      # add additive parent
      
      sign.eff <- sign(Qi[, 1])
      test.sign <- function(x){ if((x == 0)){NA} else if(x < 0){2} else{3}}
      
      par.add.id <- unlist(lapply(sign.eff, FUN = test.sign))
      Add.parent <- diag(mppData$par.per.cross[1:mppData$n.cr, par.add.id])
      
      Qi <- data.frame(Qi, Add.parent, stringsAsFactors = FALSE)
      Qi[, 1] <- abs(Qi[, 1])
      
      # add parent genotype (if given)
      
      if(!is.null(mppData$geno.par)){
        
        Par.all <- mppData$geno.par[mppData$geno.par[, 1] == QTL.list[i],
                                    5:dim(mppData$geno.par)[2]]
        Par.all <- unlist(Par.all)[Add.parent]
        
        Qi <- data.frame(Qi, Par.all, stringsAsFactors = FALSE)
        
      }
      
      # add column names
      
      if(VCOV == "h.err"){
        col.names <- c("Effect", "Std.Err", "t-test", "p-value")
      } else {col.names <- c("Effect", "Std.Err", "W-stat", "p-value")}
      
      colnames(Qi)[1:4] <- col.names
      
      Qeff.mat[[i]] <- Qi
      
    }
    
    return(Qeff.mat)
    
  } else if (Q.eff == "par"){
    
    # fill the reference parent value
    
    ref.mat <- matrix(rep(c(0, 0, 0, 1), (mppData$n.par*n.QTL)),
                      nrow = (mppData$n.par*n.QTL), byrow = TRUE)
    
    # ref names
    
    ref.names <- c()
    
    for(i in 1:length(allele_order)) {
      
      ref.names <- c(ref.names, paste0("Q", i , allele_order[[i]]))
      
    }
    
    index <- match(rownames(results), ref.names)
    ref.mat[index, ] <- results
    
    # add sign stars
    
    Sign <- sapply(ref.mat[, 4], FUN = sign.star)
    
    results <- data.frame(ref.mat, Sign, stringsAsFactors = FALSE)
    
    # add connected part references
    
    results <- data.frame(results, Con.part = unlist(con.ind),
                          stringsAsFactors = FALSE)
    
    Q.ind <- rep(paste0("Q", 1:n.QTL), each = mppData$n.par)
    Q.id <- unique(Q.ind)
    Qeff.mat <- vector(mode = "list", n.QTL)
    
    for(i in 1:n.QTL){
      
      # subset QTL
      
      Qi <- results[Q.ind == Q.id[i], ]
      rownames(Qi) <- allele_order[[i]]
      
      # add parent genotype (if given)
      
      if(!is.null(mppData$geno.par)){
        
        Par.all <- mppData$geno.par[mppData$geno.par[, 1] == QTL.list[i],
                                    5:dim(mppData$geno.par)[2]]
        
        Par.all <- unlist(Par.all)[allele_order[[i]]]
        
        Qi <- data.frame(Qi, Par.all, stringsAsFactors = FALSE)
        
      }
      
      # add column names
      
      if(VCOV == "h.err"){
        col.names <- c("Effect", "Std.Err", "t-test", "p-value")
      } else {col.names <- c("Effect", "Std.Err", "W-stat", "p-value")}
      
      colnames(Qi)[1:4] <- col.names
      
      Qeff.mat[[i]] <- Qi
      
    }
    
    return(Qeff.mat)
    
    
  } else if (Q.eff == "anc") {
    
    # form reference matrix
    
    n.allele <- lapply(X = Q.list, function(x) dim(x)[2])
    Q.ind <- rep(paste0("Q", 1:n.QTL), n.allele)
    n.allele <- sum(unlist(n.allele))
    
    ref.mat <- matrix(rep(c(0, 0, 0, 1), n.allele), nrow = n.allele,
                      byrow = TRUE)
    
    # make reference names
    
    ref.names <- c()
    
    for(i in 1:length(allele_order)) {
      
      ref.names <- c(ref.names, paste0("Q", i , allele_order[[i]]))
      
    }
    
    index <- match(rownames(results), ref.names)
    ref.mat[index, ] <- results
    rownames(ref.mat) <- ref.names
    
    # add sign stars
    
    Sign <- sapply(ref.mat[, 4], FUN = sign.star)
    
    results <- data.frame(ref.mat, Sign, stringsAsFactors = FALSE)
    
    # add connected part references
    
    results <- data.frame(results, Con.part = unlist(con.ind),
                          stringsAsFactors = FALSE)
    
    # project into parents
    
    Qeff.mat <- vector(mode = "list", n.QTL)
    ref.Q <- unique(Q.ind)
    
    for(i in 1:n.QTL){
      
      Q.mat <- results[Q.ind == ref.Q[i], ]
      A.allele <- factor(mppData$par.clu[QTL.list[i], ])
      A <- model.matrix(~ A.allele - 1)
      
      # modify column order of A
      
      x <- unlist(lapply(X = gregexpr(pattern = "A", rownames(Q.mat)),
                         FUN = function(x) x[1]))
      all.ord <- substr(rownames(Q.mat), x, nchar(rownames(Q.mat)))
      A <- A[, all.ord]
      
      # separate numeric results from con part information
      
      num.res <- as.matrix(Q.mat[, 1:4])
      num.res <- cbind(num.res, 1:dim(num.res)[1])
      
      num.res[is.na(num.res)] <- 9999
      proj.num.res <- A %*% num.res
      proj.num.res[proj.num.res == 9999] <- NA
      rownames(proj.num.res) <- mppData$parents
      
      # add the connected part information
      
      con.part.vec <- Q.mat[proj.num.res[, 5], c(5, 6), drop = FALSE]
      
      res <- data.frame(proj.num.res, con.part.vec, stringsAsFactors = FALSE)
      res <- res[order(res[, 5]), ]
      res <- res[, -5]
      
      # add parental score if provided
      
      if(!is.null(mppData$geno.par)){
        
        Par.all <- mppData$geno.par[mppData$geno.par[, 1] == QTL.list[i],
                                    5:dim(mppData$geno.par)[2]]
        
        Par.all <- unlist(Par.all)[rownames(res)]
        
        res <- data.frame(res, Par.all, stringsAsFactors = FALSE)
        
      }
      
      # add columns names
      
      if(VCOV == "h.err"){
        col.names <- c("Effect", "Std.Err", "t-test", "p-value")
      } else {col.names <- c("Effect", "Std.Err", "W-stat", "p-value")}
      
      colnames(res)[1:4] <- col.names
      
      Qeff.mat[[i]] <- res 
      
    }
    
    return(Qeff.mat)
    
    
  } else if (Q.eff == "biall"){
    
    ref.mat <- matrix(rep(c(0, 0, 0, 1), n.QTL), nrow = n.QTL, byrow = TRUE)
    ref.names <- paste0("Q", 1:n.QTL)
    index <- match(rownames(results), ref.names)
    ref.mat[index, ] <- results
    rownames(ref.mat) <- ref.names
    
    Qeff.mat <- vector(mode = "list", n.QTL)
    
    for(i in 1:n.QTL){
      
      # subset QTL
      
      Qi <- ref.mat[i, , drop = FALSE]
      
      # project into parents and add genotype score if possible
      
      if(!is.null(mppData$geno.par)){
        
        ref.mat2 <- matrix(rep(c(0, 0, 0, 1), mppData$n.par),
                           nrow = mppData$n.par, byrow = TRUE)
        
        Par.all <- mppData$geno.par[mppData$geno.par[, 1] == QTL.list[i],
                                    5:dim(mppData$geno.par)[2]]
        Par.all <- unlist(Par.all)
        
        index <- which(mppData$geno.par[, 1] == QTL.list[i])
        ref.all <- c(mppData$allele.ref[1, index, drop = FALSE])
        het.sc <- mppData$allele.ref[c(3, 4), index]
        
        ind.na <- which(is.na(Par.all))
        ind.ref <- which(Par.all == ref.all)
        ind.het <- which(((Par.all == het.sc[1])|(Par.all == het.sc[2])))
        
        ref.mat2[ind.ref, ] <- matrix(rep(Qi, length(ind.ref)),
                                      nrow = length(ind.ref), byrow = TRUE)
        
        ref.mat2[ind.ref, 1] <- 2 * ref.mat2[ind.ref, 1] 
        
        ref.mat2[ind.na, ] <- NA
        
        ref.mat2[ind.het, ] <- matrix(rep(Qi, length(ind.het)),
                                       nrow = length(ind.het), byrow = TRUE)
        
        # add the sign stars
        
        Sign <- sapply(ref.mat2[, 4], FUN = sign.star)
        ref.mat2 <- data.frame(ref.mat2, Sign, stringsAsFactors = FALSE)
        
        # add parents scores
        
        Qi <- data.frame(ref.mat2, Par.all, stringsAsFactors = FALSE)
        
        
      } else {
        
        # add Sign scores (stars)
        
        Sign <- sapply(Qi[, 4], FUN = sign.star)
        Qi <- data.frame(Qi, Sign, stringsAsFactors = FALSE)
        
      }
      
      # add columns names
      
      if(VCOV == "h.err"){
        col.names <- c("Effect", "Std.Err", "t-test", "p-value")
      } else {col.names <- c("Effect", "Std.Err", "W-stat", "p-value")}
      
      colnames(Qi)[1:4] <- col.names
      
      Qeff.mat[[i]] <- Qi
      
    }
    
    return(Qeff.mat)
    
  }
  
}