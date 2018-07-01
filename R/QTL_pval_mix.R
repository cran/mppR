################
# QTL_pval_mix #
################

# function to get the QTL decomposed genetic effect for the cross-specific
# parental and ancestral mixed models.


# QTL_pval_mix <- function(model, Q.eff, QTL.el, x, ref.name, par.names, par.clu,
#                          fct) {
#   
#   if(fct == "SIM"){
#     start.ind <- 2; end.ind <- 1
#   } else if(fct == "CIM") {
#     start.ind <- 3; end.ind <- 2}
#   
#   sign <- sign(rev(model$coefficients$fixed[1:QTL.el]))
#   pval <- asreml::wald(model)[start.ind:(QTL.el + end.ind), 4]
#   pval <- pval * sign
#   pval[pval == 0] <- 1
#   names(pval) <- ref.name
#   
#   if(Q.eff == "par"){
#     
#     pval <- pval[par.names]
#     
#   } else if (Q.eff == "anc") {
#     
#     # project into parents
#     
#     ref.all <- paste0("A.allele", par.clu[x, ])
#     pval <- pval[ref.all]
#     
#     
#   }
#   
#   return(pval)
#   
#   
# }