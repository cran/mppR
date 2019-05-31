#################
# Test mpp_proc #
#################

context("Test mpp_proc")

library(testthat)
library(mppR)

# references values (taken on the USNAM example - 24/05/2019)

# Reference QTL effects

QTL_cr <- c(2.6841679, 5.5525152, 0.6753642, 2.3775748, 4.0971491)

QTL_par <- c(0.0000000, -5.5525152, -2.6841679, -0.6753642, -4.0971491,
             -2.3775748)

QTL_anc <- c(0.000000, -1.626181, -1.626181, -3.201762, -3.201762,
             -5.552515)

QTL_biall <- c(5.207004, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000)

Ref_QTL <- list(QTL_cr, QTL_par, QTL_anc, QTL_biall)

data(mppData)

Qeff <- c('cr', 'par', 'anc', 'biall')

for(i in 1:4){
  
  res <- mpp_proc(pop.name = "USNAM", trait.name = "ULA",
                  mppData = mppData, Q.eff = Qeff[i],
                  plot.gen.eff = FALSE, CI = FALSE,
                  output.loc = tempdir(), verbose = FALSE)
  
  test_i <- paste('Q.eff =', Qeff[i])
  
  test_that(test_i, {
    
    expect_equal(object = res$QTL.effects$Qeff$Q1$Effect,
                 expected = Ref_QTL[[i]], tolerance = .0001)
    
  })
  
  
}