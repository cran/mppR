## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mppR)

## ----geno_data----------------------------------------------------------------
data("USNAM_geno", package = "mppR")
dim(USNAM_geno)
rownames(USNAM_geno)[1:6]
table(substr(rownames(USNAM_geno)[-c(1:6)], 1, 4))

## ----off_par_geno_data--------------------------------------------------------
geno.off <- USNAM_geno[7:506, ]
geno.par <- USNAM_geno[1:6, ]

## ----map_data-----------------------------------------------------------------
data("USNAM_map", package = "mppR")
head(USNAM_map)
map <- USNAM_map

## ----pheno_data---------------------------------------------------------------
data("USNAM_pheno", package = "mppR")
head(USNAM_pheno)
pheno <-  USNAM_pheno
cross.ind <- substr(rownames(pheno), 1, 4)

## ----par_per_cross------------------------------------------------------------
par.per.cross <- cbind(unique(cross.ind), rep("B73", 5),
                       rownames(geno.par)[2:6])
par.per.cross

## ----design_connectivity, fig.height = 4, fig.width = 6-----------------------
ppc_ex <- cbind(paste0("c", 1:7),
                c("PA", "PA", "PB", "PA", "PE", "PE", "PG"),
                c("PB", "PC", "PC", "PD", "PF", "PG", "PF"))

design_connectivity(ppc_ex)

## ----create.mppData-----------------------------------------------------------
mppData <- create.mppData(geno.off = geno.off, geno.par = geno.par,
                          map = map, pheno = pheno, cross.ind = cross.ind,
                          par.per.cross = par.per.cross)

## ----QC.mppData---------------------------------------------------------------
mppData <- QC.mppData(mppData = mppData, n.lim = 15, MAF.pop.lim = 0.05,
                      MAF.cr.miss = TRUE, mk.miss = 0.1,
                      gen.miss = 0.25, verbose = TRUE)

## ----IBS.mppData--------------------------------------------------------------
mppData <- IBS.mppData(mppData = mppData)


## ----IBD.mppData--------------------------------------------------------------
mppData <- IBD.mppData(mppData = mppData, het.miss.par = TRUE, type = 'RIL',
                       type.mating = 'selfing')

## ----parent_clustering2-------------------------------------------------------
data("par_clu", package = "mppR")

mppData <- parent_cluster.mppData(mppData = mppData, par.clu  = par_clu)


## ----parent_clustering, eval=FALSE, echo=TRUE---------------------------------
#  library("clusthaplo")
#  set.seed(68769)
#  
#  mppData <- parent_cluster.mppData(mppData = mppData, method = "clusthaplo",
#                                    K = 10, window = 25, plot = FALSE)
#  
#  mppData$n.anc

## ----summary_mppData----------------------------------------------------------
summary(mppData)

## ----subset_mppData-----------------------------------------------------------
mppData_sub <- subset(x = mppData, mk.list = mppData$map[, 2] == 1,
                      gen.list = sample(mppData$geno.id, 200))

## ----mpp_proc-----------------------------------------------------------------
QTL_proc <- mpp_proc(pop.name = "USNAM", trait.name = "ULA", trait = "ULA",
                     mppData = mppData, Q.eff = "anc",
                     plot.gen.eff = TRUE, N.cim = 2, thre.cof = 3,
                     win.cof = 20, window = 20, thre.QTL = 3,
                     win.QTL = 20, CI = TRUE, drop = 1.5,
                     verbose = FALSE, output.loc = tempdir())

## ----MQE----------------------------------------------------------------------
MQE <- MQE_proc(pop.name = "USNAM", trait.name = "ULA", mppData = mppData,
                Q.eff = c("par", "anc", "biall"), window = 20,
                plot.MQE = TRUE, verbose = FALSE, output.loc = tempdir())

## ----QTL_effect---------------------------------------------------------------
SIM <- mpp_SIM(mppData = mppData, Q.eff = "anc")
cofactors <- QTL_select(Qprof = SIM)
CIM <- mpp_CIM(mppData = mppData, Q.eff = "anc", cofactors = cofactors,
               plot.gen.eff = TRUE)
QTL <- QTL_select(Qprof = CIM)
gen.eff <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "anc")

summary(gen.eff, QTL = 1)

## ----QTL_profile, fig.height = 6, fig.width = 10------------------------------
plot(x = CIM, QTL = QTL, type = "l")

## ----gen_eff_plot, fig.height = 6, fig.width = 10-----------------------------
plot(x = CIM, gen.eff = TRUE, mppData = mppData, QTL = QTL, Q.eff = "anc")

## ----CV_proc, fig.height = 6, fig.width = 10----------------------------------
set.seed(89341)

CV <- mpp_CV(pop.name = "USNAM", trait.name = "ULA", mppData = mppData,
             Q.eff = "cr", her = 0.4, Rep = 2, k = 5, verbose = FALSE,
             output.loc = tempdir())

