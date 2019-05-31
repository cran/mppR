################
# Test mpp_CIM #
################

context("Test mpp_CIM")

# library(testthat)
# library(mppR)

# references values (taken on the USNAM example - 24/05/2019)

ref_ind <- c(8, 9, 11, 20, 24, 32, 37, 57, 68, 94)

ref_ind2 <- c(43, 74)

# reference -log10(p-valu)

logp_cr <- c(1.10021694, 1.68789781, 3.59980702, 0.22957343, 3.69385531,
             0.02685481, 1.43415436, 3.66469958, 1.16856072, 0.52939162)

logp_par <- c(1.10021694, 1.68789781, 3.59980702, 0.22957343, 3.69385531,
              0.02685481, 1.43415436, 3.66469958, 1.16856072, 0.52939162)

logp_anc <- c(1.3194567, 2.3176857, 4.2831107, 0.4738949, 3.6928160,
              0.1501479, 1.8364183, 4.0360087, 1.4299513, 0.6208097)

logp_biall <- c(1.90194112, 2.00673871, 2.30011066, 1.71703106, 2.19542051,
                0.40400271, 1.12271256, 0.21386840, 0.43960028, 0.01652987)

Ref_logp <- cbind(logp_cr, logp_par, logp_anc, logp_biall)

# reference QTL allele p-val

pval_cr <- c(0.05659636, 0.66696803, 0.75184962, 0.31062668, 0.63682840,
             0.26334557, -0.29658559, -0.35759213, 0.02310737, 0.04511828)

pval_par <- c(1.00000000, 1.00000000, 0.05659636, 0.66696803, 0.75184962,
              0.31062668, 0.63682840, 0.26334557, -0.29658559, -0.35759213,
              0.02310737,  0.04511828)

pval_anc <- c(1.00000000, 1.00000000, 0.09159886, 0.30683944, 0.75371835,
              0.30683944, 0.62045555, 0.91931724, 0.09159886, 0.91931724,
              0.09159886, 0.04863690)

Ref_pval <- list(pval_cr, pval_par, pval_anc)

data(mppData)

Qeff <- c('cr', 'par', 'anc', 'biall')

for(i in 1:4){
  
  CIM <- mpp_CIM(mppData = mppData, Q.eff = Qeff[i],
                 cofactors = c('L00585', 'L00194'),
                 plot.gen.eff = Qeff[i]!='biall')
  
  test_i <- paste('Q.eff =', Qeff[i])
  
  test_that(test_i, {
    
    expect_equal(object = CIM$log10pval[ref_ind], expected = Ref_logp[, i],
                 tolerance = .001)
    
  })
  
  if(Qeff[i]!='biall'){
    
    test_i <- paste0('Q_pval_', Qeff[i])
    
    test_that(test_i, {
      
      pval_i <- unlist(CIM[ref_ind2, 6:dim(CIM)[2]])
      names(pval_i) <- NULL
      
      expect_equal(object = pval_i, expected = Ref_pval[[i]], tolerance = .001)
      
    })
    
  }
  
}