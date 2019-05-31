######################
# Test mpp_back_elim #
######################

context("Test mpp_back_elim")

library(testthat)
library(mppR)

# references values (taken on the USNAM example - 24/05/2019)

# reference QTL positions entering in the model

Qcr_pos <- c(21, 35, 89)
Qpar_pos <- c(22, 53, 98)
Qanc_pos <- c(20, 36, 49)
Qbiall_pos <- c(25, 49, 63)

Q_pos_ref <- list(Qcr_pos, Qpar_pos, Qanc_pos, Qbiall_pos)

# Reference QTL staying after backward elimination.

QTL_cr <- c('L00120')

QTL_par <- c('L00410')

QTL_anc <- c('L00830')

QTL_biall <- c("L00356", "L00742")

Ref_QTL <- list(QTL_cr, QTL_par, QTL_anc, QTL_biall)

data(mppData)

Qeff <- c('cr', 'par', 'anc', 'biall')

for(i in 1:4){
  
  QTLs <- mppData$map$mk.names[Q_pos_ref[[i]]]
  
  res <- mpp_back_elim(mppData = mppData, Q.eff = Qeff[i], QTL = QTLs)
  
  test_i <- paste('Q.eff =', Qeff[i])
  
  test_that(test_i, {
    
    expect_identical(object = res$mk.names, expected = Ref_QTL[[i]])
    
  })
  
  
}