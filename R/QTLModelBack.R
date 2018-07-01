################
# QTLModelBack #
################

# Computation of the p-values of each QTL position for bacward elimination

# arguments

# mppData mppData object

# Q.list list of QTL incidence matrix

# cross.mat cross intercept QTL incidence matrix

# x model formula

# VCOV type of variance covariance structure


QTLModelBack <- function(x, mppData, trait, Q.list, cross.mat, VCOV){
  
  # 1. compute the QTL groups for asreml formula
  ##############################################
  
  if(VCOV != "h.err"){
    
    n.QTL.el <- unlist(lapply(Q.list, function(x) dim(x)[2]))
    QTL.seq <- lapply(n.QTL.el, function(x) seq(1:x))
    cum.sum <- cumsum(n.QTL.el)
    add.el <- c(0, cum.sum[-length(cum.sum)])
    QTL.seq <- mapply(function(x, y) x + y, x = QTL.seq, y = add.el,
                      SIMPLIFY = FALSE)
    
  }
  
  # 2. computation of the different models
  ########################################
  
  if(VCOV == "h.err"){
    
    an.table <- anova(lm(as.formula(x), data = Q.list))
    res <- an.table[(dim(an.table)[1] - 1), 5]
    
  } else if ((VCOV == "h.err.as") || (VCOV == "cr.err")) {
    
    # dataset <- data.frame(QTL = do.call(cbind, Q.list),
    #                       cr.mat = factor(mppData$cross.ind,
    #                                       levels = unique(mppData$cross.ind)),
    #                       trait = trait)
    # 
    # if(VCOV == "h.err.as"){ formula.R <- "~idv(units)"
    # } else if (VCOV == "cr.err") {formula.R <- "~at(cr.mat):units"}
    # 
    # model <- asreml::asreml(fixed = as.formula(x), rcov = as.formula(formula.R),
    #                 group = QTL.seq, data = dataset, trace = FALSE,
    #                 na.method.Y = "omit", na.method.X = "omit")
    # 
    # w.table <- asreml::wald(model)
    # res <- w.table[(dim(w.table)[1] - 1), 4]
    
  } else if ((VCOV == "pedigree") || (VCOV == "ped_cr.err")) {
    
    # dataset <- data.frame(QTL = do.call(cbind, Q.list),
    #                       cr.mat = factor(mppData$cross.ind,
    #                                       levels = unique(mppData$cross.ind)),
    #                       trait = trait,
    #                       genotype = mppData$geno.id)
    # 
    # if(VCOV == "pedigree"){ formula.R <- "~idv(units)"
    # } else if (VCOV == "ped_cr.err") { formula.R <- "~at(cr.mat):units"}
    # 
    # model <- asreml::asreml(fixed = as.formula(x), random = ~ ped(genotype),
    #        rcov = as.formula(formula.R), ginverse = list(genotype = ped.mat.inv),
    #        group = QTL.seq, data = dataset, trace = FALSE, na.method.Y = "omit",
    #        na.method.X = "omit")
    # 
    # w.table <- asreml::wald(model)
    # res <- w.table[(dim(w.table)[1] - 1), 4]
    
  }
  
  return(res)
  
}