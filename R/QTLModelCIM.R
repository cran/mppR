###############
# QTLModelCIM #
###############

# function to compute a single position QTL model


QTLModelCIM <- function(x, mppData, trait, cross.mat, Q.eff, VCOV,
                        cof.list, cof.part, plot.gen.eff){
  
  # 1. formation of the QTL incidence matrix
  ###########################################
  
  ### 2.2 QTL position
  
  QTL <- inc_mat_QTL(x = x, mppData = mppData, Q.eff = Q.eff, order.MAF = TRUE)
  
  QTL.el <- dim(QTL)[2] # number of QTL elements
  
  ref.name <- colnames(QTL)
  
  ### 2.1 cofactors
  
  cof.mat <- do.call(cbind, cof.list[which(cof.part[x, ])])
  
  # test if no cofactors
  
  if(is.null(cof.mat)){ cof.mat <- rep(0, length(mppData$geno.id)); cof.el <- 1
  
  } else { cof.el <- dim(cof.mat)[2] }
  
  
  # 2. model computation
  ######################
  
  ### 2.1 homogeneous residual variance error
  
  if(VCOV == "h.err"){
    
    model <- tryCatch(expr = lm(trait ~ - 1 + cross.mat + cof.mat
                                + QTL), error = function(e) NULL)
    
    if (is.null(model)){ 
      
      if(plot.gen.eff) {
        
        if(Q.eff == "cr"){ results <- c(0, rep(1, mppData$n.cr))
        
        } else { results <- c(0, rep(1, mppData$n.par)) }
        
      } else { results <- 0 }
      
    } else {
      
      
      if(!("QTL" %in% rownames(anova(model)))){ # QTL effect could not be
        # estimated probably due to
        # singularities.
        
        if(plot.gen.eff) {
          
          if(Q.eff == "cr"){ results <- c(0, rep(1, mppData$n.cr))
          
          } else { results <- c(0, rep(1, mppData$n.par)) }
          
        } else { results <- 0 }
        
      } else {
        
        results <- -log10(anova(model)$Pr[which(rownames(anova(model))=="QTL")])
        
        if(plot.gen.eff){
          
          gen.eff <- QTL_pval(mppData = mppData, model = model,
                              Q.eff = Q.eff, x = x)
          
          
          results <- c(results, gen.eff)
          
        }
        
      }
      
    }
    
    ### 2.2 HRT REML or cross-specific variance residual terms
    
  } else if ((VCOV == "h.err.as") || (VCOV == "cr.err")){
    
    # dataset <- data.frame(cof.mat = cof.mat, QTL = QTL,
    #                       cr.mat = factor(mppData$cross.ind,
    #                                       levels = unique(mppData$cross.ind)),
    #                       trait = trait)
    # 
    # colnames(dataset) <- c(paste0("cof",1:cof.el),paste0("Q",1:QTL.el),
    #                        "cr.mat", "trait")
    # 
    # formula.QTL <- paste("+", paste0("Q", 1:QTL.el), collapse = " ")
    # formula.fix <- paste("trait ~ -1 + cr.mat + grp(cof)", formula.QTL)
    # 
    # if(VCOV == "h.err.as"){ formula.R <- "~idv(units)"
    # } else if (VCOV == "cr.err") {formula.R <- "~at(cr.mat):units"} 
    # 
    # model <- tryCatch(expr = asreml::asreml(fixed = as.formula(formula.fix),
    #                                 rcov =  as.formula(formula.R),
    #                                 group = list(cof=1:cof.el),
    #                                 data=dataset, trace = FALSE,
    #                                 na.method.Y = "omit",
    #                                 na.method.X = "omit"),
    #                   error = function(e) NULL)
    
    ### 2.3 random pedigree + HVRT or + CSRT
    
  } else if ((VCOV == "pedigree") || (VCOV == "ped_cr.err")){
    
    # dataset <- data.frame(cof.mat = cof.mat, QTL = QTL,
    #                       cr.mat = factor(mppData$cross.ind,
    #                                       levels = unique(mppData$cross.ind)),
    #                       trait = trait,
    #                       genotype = mppData$geno.id)
    # 
    # colnames(dataset) <- c(paste0("cof", 1:cof.el), paste0("Q", 1:QTL.el),
    #                        "cr.mat", "trait", "genotype")
    # 
    # 
    # formula.QTL <- paste("+",paste0("Q",1:QTL.el),collapse = " ")
    # formula.fix <- paste("trait~1+grp(cof)",formula.QTL)
    # 
    # if(VCOV == "pedigree"){ formula.R <- "~idv(units)"
    # } else if (VCOV == "ped_cr.err") {formula.R <- "~at(cr.mat):units"}
    # 
    # model <- tryCatch(expr = asreml::asreml(fixed = as.formula(formula.fix),
    #                                 random = ~ ped(genotype),
    #                                 rcov = as.formula(formula.R),
    #                                 ginverse = list(genotype = ped.mat.inv),
    #                                 group = list(cof = 1:cof.el),
    #                                 data = dataset, trace = FALSE,
    #                                 na.method.Y = "omit",
    #                                 na.method.X = "omit"),
    #                   error = function(e) NULL)
    
  } 
  
  # 3. organise the results for the mixed models similar for all models
  #####################################################################
  
  # if(VCOV != "h.err"){
  #   
  #   if (is.null(model)){ 
  #     
  #     if(plot.gen.eff) {
  #       
  #       if(Q.eff == "cr"){ results <- c(0, rep(1, mppData$n.cr))
  #       
  #       } else { results <- c(0, rep(1, mppData$n.par)) }
  #       
  #     } else { results <- 0 }
  #     
  #   } else {
  #     
  #     W.stat <- sum(asreml::wald(model)[3:(QTL.el+2), 3])
  #     
  #     if(W.stat == 0){
  #       
  #       if(plot.gen.eff) {
  #         
  #         if(Q.eff == "cr"){ results <- c(0, rep(1, mppData$n.cr))
  #         
  #         } else { results <- c(0, rep(1, mppData$n.par)) }
  #         
  #       } else { results <- 0 }
  #       
  #     } else {
  #       
  #       df <- sum(asreml::wald(model)[3:(QTL.el+2), 1])
  #       
  #       pval <- pchisq(W.stat, df, lower.tail = FALSE)
  #       
  #       results <- -log10(pval)
  #       
  #       if(plot.gen.eff){
  #         
  #         gen.eff  <- QTL_pval_mix(model = model, Q.eff = Q.eff, QTL.el = QTL.el,
  #                                  x = x, ref.name = ref.name,
  #                                  par.clu = mppData$par.clu,
  #                                  par.names = mppData$parents, fct = "CIM")
  #         
  #         results  <- c(results, gen.eff)
  #         
  #       }
  #       
  #     }
  #     
  #   }
  #   
  # }
  
  return(results)
  
}