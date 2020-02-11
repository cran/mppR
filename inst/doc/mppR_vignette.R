### R code from vignette source 'mppR_vignette.rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: package_load
###################################################
library(mppR)


###################################################
### code chunk number 2: genotype_data
###################################################
data("USNAM_geno", package = "mppR")
dim(USNAM_geno)
rownames(USNAM_geno)[1:6]
table(substr(rownames(USNAM_geno)[-c(1:6)], 1, 4))


###################################################
### code chunk number 3: off_par_genotype_data
###################################################
geno.off <- USNAM_geno[7:506, ]
geno.par <- USNAM_geno[1:6, ]


###################################################
### code chunk number 4: map_data
###################################################
data("USNAM_map", package = "mppR")
head(USNAM_map)
map <- USNAM_map


###################################################
### code chunk number 5: phenotype_data
###################################################
data("USNAM_pheno", package = "mppR")
head(USNAM_pheno)
pheno <-  USNAM_pheno
cross.ind <- substr(rownames(pheno), 1, 4)


###################################################
### code chunk number 6: par.per.cross
###################################################
par.per.cross <- cbind(unique(cross.ind), rep("B73", 5),
                       rownames(geno.par)[2:6])


###################################################
### code chunk number 7: design_connectivity
###################################################
ppc_ex <- cbind(paste0("c", 1:7),
                c("PA", "PA", "PB", "PA", "PE", "PE", "PG"),
                c("PB", "PC", "PC", "PD", "PF", "PG", "PF"))

design_connectivity(ppc_ex)


###################################################
### code chunk number 8: create.mppData
###################################################
mppData <- create.mppData(geno.off = geno.off, geno.par = geno.par,
                          map = map, pheno = pheno, cross.ind = cross.ind,
                          par.per.cross = par.per.cross)


###################################################
### code chunk number 9: QC.mppData
###################################################
mppData <- QC.mppData(mppData = mppData, n.lim = 15, MAF.pop.lim = 0.05,
                      MAF.cr.miss = TRUE, mk.miss = 0.1,
                      gen.miss = 0.25, verbose = TRUE)


###################################################
### code chunk number 10: IBS.mppData
###################################################
mppData <- IBS.mppData(mppData = mppData)



###################################################
### code chunk number 11: IBD.mppData
###################################################
mppData <- IBD.mppData(mppData = mppData, het.miss.par = TRUE, type = 'RIL',
                       type.mating = 'selfing')


###################################################
### code chunk number 12: parent_clustering2
###################################################
data("par_clu", package = "mppR")

mppData <- parent_cluster.mppData(mppData = mppData, par.clu  = par_clu)



###################################################
### code chunk number 13: summary_mppData
###################################################
summary(mppData)


###################################################
### code chunk number 14: subset_mppData
###################################################
mppData_sub <- subset(x = mppData, mk.list = mppData$map[, 2] == 1,
                      gen.list = sample(mppData$geno.id, 200))


###################################################
### code chunk number 15: mpp_proc
###################################################
QTL_proc <- mpp_proc(pop.name = "USNAM", trait.name = "ULA", trait = "ULA",
                     mppData = mppData, Q.eff = "anc",
                     plot.gen.eff = TRUE, N.cim = 1, thre.cof = 3,
                     win.cof = 20, window = 20, thre.QTL = 3,
                     win.QTL = 20, CI = TRUE, drop = 1.5,
                     verbose = FALSE, output.loc = tempdir())


###################################################
### code chunk number 16: MQE (eval = FALSE)
###################################################
## MQE <- MQE_proc(pop.name = "USNAM", trait.name = "ULA", mppData = mppData,
##                 Q.eff = c("par", "anc", "biall"), window = 20, verbose = FALSE,
##                 output.loc = tempdir())


###################################################
### code chunk number 17: QTL_effect
###################################################
SIM <- mpp_SIM(mppData = mppData, Q.eff = "anc")
cofactors <- QTL_select(Qprof = SIM)
CIM <- mpp_CIM(mppData = mppData, Q.eff = "anc", cofactors = cofactors,
               plot.gen.eff = TRUE)
QTL <- QTL_select(Qprof = CIM)
gen.eff <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "anc")

summary(gen.eff, QTL = 1)


###################################################
### code chunk number 18: QTL_profile
###################################################
plot(x = CIM, QTL = QTL, type = "l")


###################################################
### code chunk number 19: gen_eff_plot
###################################################
plot(x = CIM, gen.eff = TRUE, mppData = mppData, QTL = QTL, Q.eff = "anc", main = 'QTL genetic effect plot')


###################################################
### code chunk number 20: CV_proc
###################################################
set.seed(89341)

CV <- mpp_CV(pop.name = "USNAM", trait.name = "ULA", mppData = mppData,
             Q.eff = "cr", her = 0.4, Rep = 1, k = 3, verbose = FALSE,
             output.loc = tempdir())


