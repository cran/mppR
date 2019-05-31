########################
# Test QTL_gen_effects #
########################

context("Test QTL_gen_effects")

library(testthat)
library(mppR)

# references values (taken on the USNAM example - 24/05/2019)

# Reference QTL effects.

QTL_cr <- c(2.6547645, 4.1944733, 2.2521606, 0.3315769, 1.5143095)

QTL_par <- c(0.0000000, -1.5143095, 0.3315769, -4.1944733, 2.2521606,
             -2.6547645)

QTL_anc <- c(0.000000, -2.902261, -2.902261, 1.349841, 1.349841, -2.593399)

QTL_biall <- c(0.000000, -5.658275, 0.000000, 0.000000, 0.000000, 0.000000)

Ref_QTL <- list(QTL_cr, QTL_par, QTL_anc, QTL_biall)

data(mppData)

Qeff <- c('cr', 'par', 'anc', 'biall')

for(i in 1:4){
  
  QTLs <- mppData$map$mk.names[c(25, 49, 63)]
  
  res <- QTL_gen_effects(mppData = mppData, Q.eff = Qeff[i], QTL = QTLs)
  
  test_i <- paste('Q.eff =', Qeff[i])
  
  test_that(test_i, {
    
    expect_equal(object = res$Qeff$Q1$Effect, expected = Ref_QTL[[i]],
                 tolerance = .0001)
    
  })
  
  
}