#################
# Test mpp_proc #
#################

context("Test MQE_proc")

library(testthat)
library(mppR)

# references values (taken on the USNAM example - 24/05/2019)

Ref_res <- c(0.000000, -1.732895, -1.732895, -3.064966, -3.064966, -6.127131)

data(mppData)

res <- MQE_proc(pop.name = "USNAM", trait.name = "ULA",
                mppData = mppData, Q.eff = c('cr', 'par', 'anc', 'biall'),
                output.loc = tempdir(), verbose = FALSE)
  
  test_that('MQE_test', {
    
    expect_equal(object = res$QTL.effects$Q1$Effect,
                 expected = Ref_res, tolerance = .0001)
    
  })
