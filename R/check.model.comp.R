####################
# check.model.comp #
####################

# function to check the format of the data provided and the argument consistency
# before QTL model computation (perm, SIM, CIM, QTLEffects, etc.)

# parameters

# mppData object of class mppData.

# trait numeric or character indicator to specify the trait

# Q.eff Character expression specifying the type of QTL effect.

# VCOV Character specify the variance covariance structure.

# plot.gen.eff Logical value specifying if p-value of the single QTL effect
# must be stored.

# n.cores : number of cluster

# cofactors list of cofactors

# QTL list of QTLs

# ref.par character expression indicating a reference parent.

# sum_zero Logical expression indicating if the model should be computed using
# a sum to zero constraint.

# mppData.ts mppData object of the training set

# fct specify which type of function to allow specific tests.


check.model.comp <- function(mppData = NULL, trait, Q.eff, VCOV,
                             plot.gen.eff = FALSE, n.cores = 1,
                             cofactors = NULL, QTL = NULL, ref.par = NULL,
                             sum_zero = NULL, mppData.ts = NULL,
                             fct = "XXX"){
  
  # 1. check mppData format
  #########################
  
  if(fct != "R2_pred"){
    
    if(is.null(mppData)){
      
      stop("'mppData' is not provided")
      
    } else {
      
      check_mppData(mppData = mppData, Q.eff = Q.eff)
      
    }
    
  }
  
  # 2. check trait
  ################
  
  if(fct != 'R2_pred'){ check_trait(trait = trait, mppData = mppData)
    
  } else { check_trait(trait = trait, mppData = mppData.ts)}
  
  # 3. check Q.eff argument
  #########################
  
  if (!(Q.eff %in% c("cr", "par", "anc", "biall"))){
    
    stop("'Q.eff' must be ", dQuote("cr"), ', ', dQuote("par"), ', ',
         dQuote("anc"), ' or ', dQuote("biall"))
    
  }
  
  # 4. check the VCOV argument
  ############################
  
  
  
  # test if the asreml function is present for the compuation of the mixed models
  
  if (fct != "R2"){
    
    if (!(VCOV %in% c("h.err", "h.err.as", "cr.err", "pedigree", "ped_cr.err"))){
      
      stop("'VCOV' must be ", dQuote("h.err"), ', ', dQuote("h.err.as"), ', ',
           dQuote("cr.err"), ', ', dQuote("pedigree"), ' or ',
           dQuote("ped_cr.err"))
      
    }
    
    
  } 
  
  # Warning if the user want to use mixed model for permutation
  
  if((fct=="perm") && (VCOV!="h.err")){
    
    message(paste("MESSAGE:",
                  "The determination of a significance threshold by permutation",
                  "using mixed models is technically possible but can take",
                  "a hugh amount of time! Due to limited computer power, we",
                  "advice to only use the linear model (VCOV = 'h.err')."))
    
    # ask user if he/she wants to continue
    
    text <- paste("Press [enter] to continue. If after you still want to stop",
                  "the process, Press Esc.")
    
    readline(prompt = text)
    
  }
  
  
  # 5. Consistency for parallelization.
  ####################################
  
  
  if ((n.cores > 1) && (VCOV != "h.err")){
    
    stop("parallelization is only possible for 'VCOV' = ", dQuote("h.err")) 
    
    
  }
  
  
  # 6. other checks according to specific type of functions
  #########################################################
  
  if((fct == "SIM")||(fct == "CIM")) {
    
    ### 5.1 test that if the user wants to fit a bi-allelic model, the
    # plot.gen.eff is not activated
    
    if((Q.eff == "biall") && plot.gen.eff) {
      
      stop("the estimation using 'plot.gen.eff' = TRUE is not available for the bi-allelic model")
      
    }
    
    if(fct == "CIM"){
      
      if(is.null(cofactors)) {
        
        stop("'cofactors' is not provided")
        
      }
      
    }
    
  }
  
  
  
  ### test the format of the QTL list introduce for backward elimination, R2 and
  # genetic effects estimation
  
  if(fct %in% c("back", "R2", "R2_pred", "QTLeffects")){
    
    if(is.null(QTL)){
      
      stop("'QTL' does not contain any QTL position")
      
    }
    
    if(is.character(QTL)) {
      
      if(sum(!(QTL %in% mppData$map[, 1])) != 0){
        
        wrong.QTL <- QTL[!(QTL %in% mppData$map[, 1])]
        
        pbQTL <- paste(wrong.QTL, collapse = ", ")
        message <- sprintf(ngettext(length(wrong.QTL),
                                    "'QTL' position %s is not present in 'mppData'",
                                    "'QTL' positions %s are not present in 'mppData'"),
                           pbQTL)
        
        stop(message)
        
      }
      
    } else { # the list of QTL is not character test if QTLlist format
      
      stopifnot(inherits(QTL, "QTLlist"))
      
    }
    
    
    if (fct == "QTLeffects"){
    
      ### Test the compatibility of the reference parents
        
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
      
      ### Test correct configuration if sum_zero = TRUE
      
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
  
}