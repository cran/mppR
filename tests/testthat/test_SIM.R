################
# Test mpp_SIM #
################

context("Test mpp_SIM")

library(testthat)
library(mppR)

# references values (taken on the USNAM example - 24/05/2019)

ref_ind <- c(6, 7, 13, 19, 22, 32, 35, 54, 71, 93)

ref_ind2 <- c(24, 64)

# reference -log10(p-valu)

logp_cr <- c(1.2122360, 0.6667540, 5.8191649, 1.4675994, 3.3868153,
            0.2748310, 1.0899923, 0.4067160, 0.3809673, 0.5170291)

logp_par <- c(1.2122360, 0.6667540, 5.8191649, 1.4675994, 3.3868153,
             0.2748310, 1.0899923, 0.4067160, 0.3809673, 0.5170291)

logp_anc <- c(0.7431981, 0.4043566, 6.7205350, 1.4675994, 3.3868153,
             0.2748310, 1.0899923, 0.4703100, 0.3563172, 0.6406431)

logp_biall <- c(0.3871609, 0.8074712, 3.1787436, 0.9081235, 1.7338572,
               0.1745960, 0.5331881, 1.5266279, 0.5902561, 0.5324062)

Ref_logp <- cbind(logp_cr, logp_par, logp_anc, logp_biall)

# reference QTL allele p-val

pval_cr <- c(-0.0430084110, -0.3288853113, -0.0003377911, 0.0360343477,
             0.0652087216, -0.4123434221, -0.6775865568, -0.0001472831,
             -0.1253379854, 0.9016560046)

pval_par <- c(1.0000000000, 1.0000000000, -0.0430084110, -0.3288853113,
              -0.0003377911, 0.0360343477, 0.0652087216, -0.4123434221,
              -0.6775865568, -0.0001472831, -0.1253379854, 0.9016560046)

pval_anc <- c(1.0000000000, 1.0000000000, -0.0430084110, -0.5418191542,
              -0.0003377911, 0.0359573213, 0.0652087216, -0.4121505542,
              -0.6775865568, -0.0001463367, -0.1253379854, -0.5418191542)

Ref_pval <- list(pval_cr, pval_par, pval_anc)

data(mppData)

Qeff <- c('cr', 'par', 'anc', 'biall')

for(i in 1:4){
  
  SIM <- mpp_SIM(mppData = mppData, Q.eff = Qeff[i],
                 plot.gen.eff = Qeff[i]!='biall')
  
  test_i <- paste('Q.eff =', Qeff[i])
  
  test_that(test_i, {
    
    expect_equal(object = SIM$log10pval[ref_ind], expected = Ref_logp[, i],
                 tolerance = .001)
    
  })
  
  if(Qeff[i]!='biall'){
    
    test_i <- paste0('Q_pval_', Qeff[i])
    
    test_that(test_i, {
      
      pval_i <- unlist(SIM[ref_ind2, 6:dim(SIM)[2]])
      names(pval_i) <- NULL
      
      expect_equal(object = pval_i, expected = Ref_pval[[i]], tolerance = .001)
      
    })
    
  }
  
}