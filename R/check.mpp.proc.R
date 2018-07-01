##################
# check.mpp.proc #
##################

# function to check the format of all elements introduced in the function
# mpp.proc

check.mpp.proc <- function(mppData, trait, Q.eff, VCOV, plot.gen.eff = FALSE,
                           ref.par = NULL, sum_zero = NULL, n.cores = 1,
                           output.loc){
  
  # 1. test the validity of the provided path to store the results
  
  if(!file.exists(output.loc)){
    
    stop("'output.loc' is not a valid path")
    
  }
  
  # 2. check mppData format
  
  check_mppData(mppData = mppData, Q.eff = Q.eff)
  
  # 3. check trait format
  
  check_trait(trait = trait, mppData = mppData)
  
  # 4. check Q.eff argument
  
  
  if (!(Q.eff %in% c("cr", "par", "anc", "biall"))){
    
    stop("'Q.eff' must be ", dQuote("cr"), ', ', dQuote("par"), ', ',
         dQuote("anc"), ' or ', dQuote("biall"))
    
  }
  
  # 5. check the VCOV argument
  
  
  if (!(VCOV %in% c("h.err", "h.err.as", "cr.err", "pedigree", "ped_cr.err"))){
    
    stop("'VCOV' must be ", dQuote("h.err"), ', ', dQuote("h.err.as"), ', ',
         dQuote("cr.err"), ', ', dQuote("pedigree"), ' or ',
         dQuote("ped_cr.err"))
    
  }
  
  # 6. test if the asreml function is present for the compuation of the
  # mixed models
  
    
    # if (VCOV != "h.err"){
    #   
    #   test <- requireNamespace(package = 'asreml', quietly = TRUE)
    #   
    #   if(!test){
    #     
    #     stop(paste("To use this type of VCOV, you must have access to the asreml",
    #                "function from the asreml-R package."))
    #     
    #   }
    #   
    # }

    
  
  # 7. Consistency for parallelization.
  
  
  if ((n.cores > 1) && (VCOV != "h.err")){
    
    stop("parallelization is only possible for 'VCOV' = ", dQuote("h.err"))  
    
    
  }
  
  
    # 8. Test that if the user wants to fit a bi-allelic model, the
    # plot.gen.eff is not activated
    
    if((Q.eff == "biall") && plot.gen.eff) {
      
      stop("the estimation using 'plot.gen.eff' = TRUE is not available for the bi-allelic model")
      
      
    }
  
  # 9. Check the argument ref.par
  
  if(!is.null(ref.par)){
    
    # test that there is only one reference parent and one connected part.
    
    if(length(ref.par) !=1){
      
      stop("'ref.par' must be of length 1")
      
    }
    
    nb.con.part <- length(design_connectivity(mppData$par.per.cross,
                                              plot_des = FALSE))
    
    if(nb.con.part > 1){
      
      stop("'ref.par' can only be used if the MPP has a unique connected part")
      
    }
    
    # test that reference parent is present in the list of parents.
    
    if(!(ref.par %in% mppData$parents)){
      
      par_list <- paste(mppData$parents, collapse = ", ")
      
      stop("'ref.par' must be one of: ", par_list)
      
    }
    
  }
  
  # 10. Check the argument sum_zero
  
  if(!is.null(sum_zero)){
    
    if(sum_zero){
      
      if(!(Q.eff %in% c('par', 'anc'))){
        
        stop("you can only use the sum to zero constraint for the parental ",
             "or the ancestral model")
        
      }
      
      if(VCOV != 'h.err'){
        
        stop("you can only use the sum to zero constraint for the ",
             "homogeneous error term model")
        
      }
      
    }
    
    
  }
    
}