# mppR 1.5.0
* Modification replacement of the functions QTL_effect_main_QxE by QTL_effect_main_QEI and QTL_effect_QxEC by QTL_effect_main_QxEC to update the procedure to estimate the significance of the QTL by environment interaction of the parental allelic effects and the estimation of the sensitivity slope
* Correction of a bug to calculate the significance of the QTLxE effect using the Wald test

# mppR 1.4.0
* Addition of function to perform multi-parent multi-environment QTL detection (mppGE_SIM, mppGE_CIM, mppGE_proc, etc.)
* Addition of a new vignette to illustrate the multi-parent multi-environment QTL detection
* Extra function (mppData_add_pheno, mppData_mdf_pheno) to manipulate mppData objects especially modify the phenotypic data

# mppR 1.3.0
* Correction of few problems to load BCsFt populations (IBD.mppData).
* Addition of new options in MQE_proc: possibility to select a reference parent
effect calculation, possibility to calculate a confidence interval for the QTL positions
* Addition of new functions (QTL_forward, mpp_forward) to realize a MPP QTL detection
using a forward regression strategy

# mppR 1.2.1
* Correction due to the fact that package synbreed was removed from CRAN
* The possibility to do marker imputation was suppressed.
* Small modification in the genetic effect plot. Inversion of color

# mppR 1.2.0
* Corrections of bugs affecting the presentation of the QTL additive effects in: QTL_gen_effects(), MQE_gen_effects(), and the QTL reports of mpp_proc() and MQE_proc()
* Review of the vignette, small modifications in the documentation and in the error messages
* Stop to update the branch 'mppR_clusthaplo'. Stay two branches: 'master' containing the largest version of the package with mixed models and direct call to clusthaplo, and 'mppR_CRAN' a reduced version without the mixed models and the call to clusthaplo. 

# mppR 1.1.11
* Put dependency to the 'synbreed' package as suggested and not depends

# mppR 1.1.10
* version of the package contained on branch 'mppR_CRAN' published on CRAN

# mppR 1.1.8
* submission to CRAN
* Revision: suppress the dependency to clusthaplo. Form two branches on Github: 'mppR_clusthaplo' (contain function to call clusthaplo) and 'mppR_CRAN'

# mppR 1.1.7
* publish mppR on Github