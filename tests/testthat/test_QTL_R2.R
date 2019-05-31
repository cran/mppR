###############
# Test QTL_R2 #
###############

context("Test QTL_R2")

library(testthat)
library(mppR)

# references values (taken on the USNAM example - 24/05/2019)

# Reference QTL partial R2.

QTL_cr <- c(4.896388, 1.948077, 2.534107)

QTL_par <- c(4.896388, 1.948077, 2.534107)

QTL_anc <- c(4.084807, 1.851781, 2.623775)

QTL_biall <- c(1.2351597, 0.8273011, 0.4149436)

Ref_QTL <- list(QTL_cr, QTL_par, QTL_anc, QTL_biall)

data(mppData)

Qeff <- c('cr', 'par', 'anc', 'biall')

for(i in 1:4){
  
  QTLs <- mppData$map$mk.names[c(25, 49, 63)]
  
  res <- QTL_R2(mppData = mppData, Q.eff = Qeff[i], QTL = QTLs)
  
  test_i <- paste('Q.eff =', Qeff[i])
  
  test_that(test_i, {
    
    part_R2_i <- res$part.R2.diff
    names(part_R2_i) <- NULL
    
    expect_equal(object = part_R2_i, expected = Ref_QTL[[i]],
                 tolerance = .0001)
    
  })
  
  
}