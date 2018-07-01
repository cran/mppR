#############
# check.MQE #
#############

# function to check the format of the data provided to MQE functions.

# parameters

# mppData:  objects of class mppData.

# Q.eff Character vector specifying the type of QTL effect.

# VCOV Character specify the variance covariance structure.

# n.cores: number of cluster

# cluster cluster to run the function in parallel


check.MQE <- function(mppData = NULL, trait, Q.eff, VCOV, cofactors = NULL,
                      cof.Qeff, n.cores = 1, QTL = NULL, output.loc,
                      fct = "XXX"){
  
  # 1. check mppData format
  #########################
  
  check_mppData(mppData = mppData, Q.eff = Q.eff)
  
  # 2. check trait
  ################
  
  check_trait(trait = trait, mppData = mppData)
  
  # 2. check Q.eff argument
  #########################
  
  if((fct == "proc") | (fct == "forward")){
    
    if (length(Q.eff) <= 1) {
      message("you should provide at least 2 type of QTL effect for 'Q.eff'")
    }
    
  }
  
  
  test.Qeff <- Q.eff %in% c("cr", "par", "anc", "biall")
  
  if(sum(!test.Qeff) != 0 ){
    
    wrong.Qeff <- Q.eff[!test.Qeff]
    pbQeff <- paste(wrong.Qeff, collapse = ", ")
    message <- sprintf(ngettext(length(wrong.Qeff),
                                "'Q.eff' element %s is not valid. Only cr, par, anc or biall are allowed",
                                "'Q.eff' elements %s are not valid. Only cr, par, anc or biall are allowed"),
                                pbQeff)
      
    stop(message)
    
  }
  
  
  # 5. check the VCOV argument
  ############################
  
  
  # test if the asreml function is present for the compuation of the mixed models
  
  if(fct != "R2") {
    
    if (!(VCOV %in% c("h.err", "h.err.as", "cr.err", "pedigree", "ped_cr.err"))){
      
      stop("'VCOV' must be ", dQuote("h.err"), ', ', dQuote("h.err.as"), ', ',
           dQuote("cr.err"), ', ', dQuote("pedigree"), ' or ',
           dQuote("ped_cr.err"))
      
    }
    
    
  }
  
  # Warning if the user want to use mixed model for
  
  if((fct == "forward") && (VCOV!="h.err")){
    
    message(paste("MESSAGE: The determination of a MQE model",
                  "using mixed models is technically possible but can take",
                  "a lot of time!"))
    
    text <- paste("Press [enter] to continue. If after you still want to stop",
                  "the process, Press Esc.")
    
    readline(prompt = text)
    
  }
  
  
  
  # 6. Consistency for parallelization.
  ####################################
  
  
  if ((n.cores > 1) && (VCOV != "h.err")){
    
    stop("parallelization is only possible for 'VCOV' = ", dQuote("h.err")) 
    
    
  }
  
  
  # X. other checks according to specific type of functions
  #########################################################
  
  
  ### test the format of the QTL list introduce for backward elimination, R2 and
  # genetic effects estimation
  
  if(fct %in% c("R2", "QTLeffects")){
    
    if(is.null(QTL)){
      
      stop("'QTL' does not contain any QTL position")
      
    }
    
    if(is.character(QTL)) {
      
      if(sum(!(QTL %in% mppData$map[, 1])) != 0){
        
        wrong.QTL <- QTL[!(QTL %in% mppData$map[, 1])]
        
        pbQTL <- paste(wrong.QTL, collapse = ", ")
        message <- sprintf(ngettext(length(wrong.QTL),
                                    "'QTL' position %s is not present in 'Qprof'",
                                    "'QTL' positions %s are not present in 'Qprof'"),
                           pbQTL)
        
        stop(message)
        
      }
      
    } else {
      
      stop("'QTL' is not character")
      
    }
    
  }
  
  if(fct == "CIM"){
    
    if(is.null(cofactors)){
      
      stop("'cofactors' is not provided")
      
    }
    
  }
  
  if(fct == "proc"){ # test the validity of the path
    
    if(!file.exists(output.loc)){
      
      stop("'output.loc' is not a valid path")
      
    }
    
  }
  
  }
