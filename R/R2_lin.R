##########
# R2_lin #
##########

# function to compute the R squared of a list of QTL

# argument

# mppData mppData

# QTL list of QTL

# adjust should the R squared value be adjusted or not

R2_lin <- function(mppData, trait, QTL){
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  # remove non complete observations
  
  dataset <- cbind(trait, cross.mat, QTL)
  index <- complete.cases(dataset)
  trait <- trait[index]
  cross.mat <- cross.mat[index, , drop = FALSE]
  QTL <- QTL[index, , drop = FALSE]
  
  
  model.full <- tryCatch(lm(trait~-1 + cross.mat + QTL),
                         error=function(e) NULL)
  
  model.red <- tryCatch(lm(trait~-1+cross.mat),
                        error=function(e) NULL)
  
  if(is.null(model.full) || is.null(model.red)){R2 <- R2.adj <- NA
  
  } else {
    
    RSS.full <-  anova(model.full)$S[3]
    RSS.red <- anova(model.red)$S[2]
    
    df.full <- anova(model.full)$Df[3]
    df.red <- anova(model.red)$Df[2]
    
    R2 <- 100 * (1-(RSS.full/RSS.red))
    
    R2.adj <- 100 * (1 - ((RSS.full/df.full)/(RSS.red/df.red)))
    
  }
  
  return(list(R2 = R2, R2.adj = R2.adj))
  
}