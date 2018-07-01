################
# QTLModelPerm #
################

# function to compute a single position QTL model

QTLModelPerm <- function(x, mppData, trait, cross.mat, Q.eff, VCOV){
  
  # 1. formation of the QTL incidence matrix
  ###########################################
  
  QTL <- inc_mat_QTL(x = x, mppData = mppData, Q.eff = Q.eff, order.MAF = TRUE)
  
  # 2. model computation
  ######################
  
  ### 2.1 homogeneous residual variance error
  
  if(VCOV == "h.err"){
    
    model <- tryCatch(expr = lm(trait ~ -1 + cross.mat + QTL),
                      error = function(e) NULL)
    
    if (is.null(model)){ results <- 0
    
    } else { results <- -log10(anova(model)$Pr[2]) }
    
    ### 2.2 HRT REML or cross-specific variance residual terms
    
  } else if ((VCOV == "h.err.as") || (VCOV == "cr.err")){
    
    # dataset <- data.frame(QTL = QTL,
    #                       cr.mat = factor(mppData$cross.ind,
    #                                       levels = unique(mppData$cross.ind)),
    #                       trait = trait)
    # 
    # if(VCOV == "h.err.as"){ formula.R <- "~idv(units)"
    # } else if (VCOV == "cr.err") {formula.R <- "~at(cr.mat):units"}
    # 
    # 
    # model <- tryCatch(expr = asreml::asreml(fixed = trait ~ -1 + cr.mat + grp(QTL),
    #                                 rcov =  as.formula(formula.R),
    #                                 group = list(QTL = 1:dim(QTL)[2]),
    #                                 data=dataset, trace = FALSE,
    #                                 na.method.Y = "omit",
    #                                 na.method.X = "omit"),
    #                   error = function(e) NULL)
    # 
    # 
    # ### 2.3 random pedigree + HVRT or + CSRT
    
  } else if ((VCOV == "pedigree") || (VCOV == "ped_cr.err")){
    
    # # compose the dataset for the asreml function
    # 
    # dataset <- data.frame(QTL = QTL, trait = trait,
    #                       cr.mat = factor(mppData$cross.ind,
    #                                       levels = unique(mppData$cross.ind)),
    #                       genotype = mppData$geno.id)
    # 
    # if(VCOV == "pedigree"){ formula.R <- "~idv(units)"
    # } else if (VCOV == "ped_cr.err") {formula.R <- "~at(cr.mat):units"}
    # 
    # model <- tryCatch(expr = asreml::asreml(fixed = trait ~ 1 + grp(QTL),
    #                                 random = ~ ped(genotype),
    #                                 rcov =  as.formula(formula.R),
    #                                 group = list(QTL=1:dim(QTL)[2]),
    #                                 ginverse = list(genotype = ped.mat.inv),
    #                                 data = dataset, trace = FALSE,
    #                                 na.method.Y = "omit",
    #                                 na.method.X = "omit"),
    #                   error = function(e) NULL)
    
  } 
  
  # if(VCOV != "h.err"){
  #   
  #   if (is.null(model)){ results <- 0
  #   
  #   } else { pval <- pchisq(asreml::wald(model)[2, 3], asreml::wald(model)[2, 1],
  #                           lower.tail = FALSE)
  #   results <- -log10(pval)
  #   
  #   }
  #   
  # }
  
  return(results)
  
}