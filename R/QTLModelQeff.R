################
# QTLModelQeff #
################

# computation of the model for the QTL genetic effects

# arguments

# mppData : data object

# Q.list: list of QTL positions

# VCOV: variance covariance structure

# trait: trait values

# cross.mat : cross intercept incidence matrix


QTLModelQeff <- function(mppData, trait, cross.mat, Q.list, VCOV){
  
  
  if(VCOV == "h.err"){
    
    n.QTL <- length(Q.list)
    
    formula <- paste("trait ~ -1 + cross.mat + ",
                     paste(paste0("Q", 1:n.QTL), collapse = "+"))
    
    model <- lm(as.formula(formula), data = Q.list)
    
    
  } else if ((VCOV == "h.err.as") || (VCOV == "cr.err")) {
    
    # n.QTL <- length(Q.list)
    # 
    # dataset <- data.frame(do.call(cbind, Q.list),
    #                       cr.mat = factor(mppData$cross.ind,
    #                                       levels = unique(mppData$cross.ind)),
    #                       trait = trait)
    # 
    # Q.names <- function(x, Q.list){
    #   paste0("Q", x, attr(Q.list[[x]], "dimnames")[[2]])
    # }
    # 
    # names.QTL <- unlist(lapply(X = 1:n.QTL, FUN = Q.names, Q.list = Q.list))
    # names.QTL <- make.names(names.QTL)
    # 
    # colnames(dataset)[1:length(names.QTL)] <- names.QTL
    # 
    # fbegin <- "trait ~ -1 + cr.mat +"
    # 
    # formula <- paste(fbegin, paste(names.QTL, collapse = "+"))
    # 
    # if(VCOV == "h.err.as"){ formula.R <- "~idv(units)"
    # } else if (VCOV == "cr.err") {formula.R <- "~at(cr.mat):units"}
    # 
    # model <- asreml::asreml(fixed = as.formula(formula), rcov = as.formula(formula.R),
    #                 data = dataset, trace = FALSE, na.method.Y = "omit",
    #                 na.method.X = "omit")
    
    
  } else if ((VCOV == "pedigree") || (VCOV == "ped_cr.err")) {
    
    # n.QTL <- length(Q.list)
    # 
    # dataset <- data.frame(do.call(cbind, Q.list),
    #                       cr.mat = factor(mppData$cross.ind,
    #                                       levels = unique(mppData$cross.ind)),
    #                       trait = trait,
    #                       genotype = mppData$geno.id)
    # 
    # Q.names <- function(x, Q.list){
    #   paste0("Q", x, attr(Q.list[[x]], "dimnames")[[2]])
    # }
    # 
    # names.QTL <- unlist(lapply(X = 1:n.QTL, FUN = Q.names, Q.list = Q.list))
    # names.QTL <- make.names(names.QTL)
    # 
    # colnames(dataset)[1:length(names.QTL)] <- names.QTL
    # 
    # fbegin <- "trait ~ 1 +"
    # 
    # formula <- paste(fbegin, paste(names.QTL, collapse = "+"))
    # 
    # 
    # if(VCOV == "pedigree"){ formula.R <- "~idv(units)"
    # } else if (VCOV == "ped_cr.err") { formula.R <- "~at(cr.mat):units"}
    # 
    # model <- asreml::asreml(fixed = as.formula(formula), random = ~ ped(genotype),
    #                 rcov = as.formula(formula.R),
    #                 ginverse = list(genotype = ped.mat.inv),
    #                 data = dataset, trace = FALSE, na.method.Y = "omit",
    #                 na.method.X = "omit")
    
    
  }
  
  return(model)
  
}