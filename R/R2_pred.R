###########
# R2_pred #
###########

# function to compute the predicted R squared in a validation set

# y.vs: observed trait value in the validation set

# B.ts: list of QTL genetic effects estimated in the training set

# Q.list: list of QTL incidence matrices

# her: single value representing the average cross heritability or vector of
# within cross heritability.


R2_pred <- function(mppData.vs, y.vs, B.ts, Q.list, her) {
  
  X.vs <- as.matrix(do.call(cbind, Q.list))
  
  # use only complete case informations
  
  dataset <- cbind(y.vs, X.vs)
  cross.ind <- mppData.vs$cross.ind
  index <- complete.cases(dataset)
  
  y.vs <- dataset[index, 1, drop = FALSE]
  X.vs <- dataset[index, 2:dim(dataset)[2], drop = FALSE]
  cross.ind <- cross.ind[index]
  n.cr <- table(factor(cross.ind, levels = unique(cross.ind)))
  
  B.ts <- unlist(B.ts)
  B.ts[is.na(B.ts)] <- 0
  
  ############# option here to detect outliers in Beta
  
  y.vs.hat <- X.vs %*% B.ts
  
  dataset <- data.frame(y.vs, y.vs.hat)
  
  
  with.cross.cor <- function(x){
    if((length(unique(x[, 1])) == 1) || (length(unique(x[, 2])) == 1)){ 0
    } else {100 * ((cor(x[, 1], x[, 2])^2))}
  }
  
  dataset.cr <- split(x = dataset,
                      f = factor(cross.ind, levels = unique(cross.ind)))
  
  R2.cr <- unlist(lapply(X = dataset.cr, FUN =  with.cross.cor))
  
  if(length(her) > 1){ her <- her[unique(mppData.vs$cross.ind)
                                  %in% unique(cross.ind)] }
  
  R2.cr <- R2.cr/her
  
  data.pvs <- cbind(R2.cr, n.cr)
  data.pvs <- data.pvs[complete.cases(data.pvs), ]
  
  R2.av <- sum(data.pvs[, 1] * data.pvs[, 2])/sum(data.pvs[, 2])
  
  return(list(R2.av, R2.cr))
  
}