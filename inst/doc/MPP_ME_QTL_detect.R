## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, fig.height = 4, fig.width = 6-------------------------------------
library(mppR)
data(mppData_GE)
design_connectivity(par_per_cross = mppData_GE$par.per.cross)

## ----SIM, fig.height = 5, fig.width = 9---------------------------------------
SIM <- mppGE_SIM(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'), ref_par = 'UH007')
plot(x = SIM)

## ----cofactors----------------------------------------------------------------
cofactors <- QTL_select(Qprof = SIM, threshold = 4, window = 50)

## ----CIM----------------------------------------------------------------------
CIM <- mppGE_CIM(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
                 VCOV = "UN", VCOV_data = "unique", cofactors = cofactors,
                 window = 20)

## ----QTLs---------------------------------------------------------------------
QTL <- QTL_select(Qprof = SIM, threshold = 4, window = 20)

## ----QTL_effects--------------------------------------------------------------
Qeff <- QTL_effect_GE(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
                       QTL = QTL)
Qeff$QTL_1

## ----QTL_effects_2------------------------------------------------------------
Qeff <- QTL_effect_GE(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
                       QTL = QTL, ref_par = 'F2')
Qeff$QTL_1

## ----QTL_R2-------------------------------------------------------------------
QR2 <- QTL_R2_GE(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'), QTL = QTL)

QR2$glb.adj.R2
QR2$part.adj.R2.diff


## ----Q_prof, fig.height = 5, fig.width = 9------------------------------------
plot(x = CIM)

## ----Q_eff_plot, fig.height = 5, fig.width = 9--------------------------------
plot_allele_eff_GE(mppData = mppData_GE, nEnv = 2,
                   EnvNames = c('CIAM', 'TUM'), Qprof = CIM,
                   QTL = QTL, text.size = 14)

## ----mppGE_proc---------------------------------------------------------------
MPP_GE_QTL <- mppGE_proc(pop.name = 'EUNAM', trait.name = 'DMY',
                         mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
                         n.cores = 1, verbose = FALSE, output.loc = tempdir())

## ----QTLxEC, eval=FALSE-------------------------------------------------------
#  EC <- matrix(c(180, 310, 240, 280), 4, 1)
#  rownames(EC) <- c('CIAM', 'TUM', 'INRA', 'KWS')
#  colnames(EC) <- 'cum_rain'
#  
#  Qeff <- QTL_effect_QxEC(mppData = mppData_GE,
#                          trait = c('DMY_CIAM', 'DMY_TUM', 'DMY_INRA_P', 'DMY_KWS'),
#                          env_id = c('CIAM', 'TUM', 'INRA', 'KWS'),
#                          QTL = QTL, EC = EC)
#  
#  Qeff$Qeff_EC$QTL1

## ----QTLxEC_plot, eval=FALSE--------------------------------------------------
#  pl <- plot_QTLxEC(Qeff, Q_id = 2, RP = "UH007", EC_id = 'cum rain', trait_id = 'DMY')
#  pl

## ----comp_res_table, echo = FALSE---------------------------------------------
library(knitr)
tab <- data.frame(Population = c('BCNAM Grinkan', 'BCNAM Kenin-Keni', 'BCNAM Lata' ,'EUNAM Dent', 'EUNAM Flint', 'Breeding pop1', 'Breeding pop2'),
                  Ngeno = c(1598, 575, 896, 841, 811, 2071, 820),
                  Nmarker = c(51545, 51545, 51545, 18621, 5949, 1812, 1760),
                  Nenv = c(4, 4, 3, 4, 6, 3, 5),
                  Ncore = c(4, 4, 4, 1, 1, 1, 1),
                  `SIM(m)` = c(20, 1.5, 3, 18, 15, 5, 9),
                  `CIM(h)` = c(11, 0.17, 0.33, 9.5, 16, 1.3, 10.2),
                  `QTL_effect(m)` = c(15, 1.5, 2, 19, 150, 9, 33),
                  `Total(h)` = c(11.6, 0.25, 0.4, 10.2, 19, 1.5, 10.8)
                  )
kable(tab, caption = 'Computation time examples',
      col.names = c('Populations', 'Ngeno', 'Nmarker', 'Nenv', 'Ncore', 'SIM [m]', 'CIM [h]', 'QTL_effects [m]', 'Total [h]'))

