###################
# check_mod_mppGE #
###################

# Function to check the argument to compute a mpp GxE model

check_mod_mppGE <- function(mppData, trait, Q.eff, VCOV, CIM=FALSE,
                            cofactors=NULL, QTL_ch=FALSE, QTL=NULL, ref_par=NULL,
                            GE=TRUE, fast = FALSE, plot_data=NULL,
                            exp_des_form=NULL){

 # 1. check mppData
 #####

  if(!is_mppData(mppData)) {

    stop("'mppData' must be of class ", dQuote("mppData"))

  }

  if((mppData$status != 'IBD') && (mppData$status != 'complete')){

    stop("'mppData' is not complete. Use first all processing ",
         "functions in the specified order: QC.mppData, IBS.mppData, ",
         "IBD.mppData, and optionally parent_cluster.mppData for the ancestral ",
         "model")

  }

  #####

  # 2. check Q.eff
  ################

  if(!is.character(Q.eff)){

    stop("'Q.eff' must be character")

  }

  if (!(Q.eff %in% c("cr", "par", "anc", "biall"))){

    stop("'Q.eff' must be ", dQuote("cr"), ', ', dQuote("par"), ', ',
         dQuote("anc"), ' or ', dQuote("biall"))

  }

  #########

  # 3. check cofactors for CIM
  ########

  if(CIM){

    if(is.null(cofactors)) {

      stop("'cofactors' is not provided")

    }

    if(!(is.character(cofactors) || inherits(cofactors, "QTLlist"))){

      stop("'cofactors' must either be a character vector or an object of class QTLlist")

    }
    
    if("QTLlist" %in% class(cofactors)){
      
      t_cof <- !(cofactors$mk.names %in% mppData$map$mk.names)
      
      if(any(t_cof)){
        
        cof_miss <- cofactors$mk.names[t_cof]
        
        m_err <- paste("The following cofactors:", paste(cof_miss, collapse = ", "),
                       "are not present in the map of the mppData object.")
        
        stop(m_err)
        
      }
      
    } else if(is.character(cofactors)){
      
      t_cof <- !(cofactors %in% mppData$map$mk.names)
      
      if(any(t_cof)){
        
        cof_miss <- cofactors[t_cof]
        
        m_err <- paste("The following cofactors:", paste(cof_miss, collapse = ", "),
                       "are not present in the map of the mppData object.")
        
        stop(m_err)
        
      }
      
    }

  }

  
  #######
  
  # 4. check QTL
  ########
  
  if(CIM){
    
    if(is.null(cofactors)) {
      
      stop("'cofactors' is not provided")
      
    }
    
    if(!(is.character(cofactors) || inherits(cofactors, "QTLlist"))){
      
      stop("'cofactors' must either be a character vector or an object of class QTLlist")
      
    }
    
    if("QTLlist" %in% class(cofactors)){
      
      t_cof <- !(cofactors$mk.names %in% mppData$map$mk.names)
      
      if(any(t_cof)){
        
        cof_miss <- cofactors$mk.names[t_cof]
        
        m_err <- paste("The following cofactors:", paste(cof_miss, collapse = ", "),
                       "are not present in the map of the mppData object.")
        
        stop(m_err)
        
      }
      
    } else if(is.character(cofactors)){
      
      t_cof <- !(cofactors %in% mppData$map$mk.names)
      
      if(any(t_cof)){
        
        cof_miss <- cofactors[t_cof]
        
        m_err <- paste("The following cofactors:", paste(cof_miss, collapse = ", "),
                       "are not present in the map of the mppData object.")
        
        stop(m_err)
        
      }
      
    }
    
  }
  
  
  #######

  # 5. check QTL for back_elim, R2 or QTL_effect
  #############

  if(QTL_ch){

    if(is.null(QTL)) {
      
      stop("'QTL' is not provided")
      
    }
    
    if(!(is.character(QTL) || inherits(QTL, "QTLlist"))){
      
      stop("'QTL' must either be a character vector or an object of class QTLlist")
      
    }
    
    if("QTLlist" %in% class(QTL)){
      
      t_cof <- !(QTL$mk.names %in% mppData$map$mk.names)
      
      if(any(t_cof)){
        
        cof_miss <- QTL$mk.names[t_cof]
        
        m_err <- paste("The following QTL:", paste(cof_miss, collapse = ", "),
                       "are not present in the map of the mppData object.")
        
        stop(m_err)
        
      }
      
    } else if(is.character(QTL)){
      
      t_cof <- !(QTL %in% mppData$map$mk.names)
      
      if(any(t_cof)){
        
        cof_miss <- QTL[t_cof]
        
        m_err <- paste("The following QTL:", paste(cof_miss, collapse = ", "),
                       "are not present in the map of the mppData object.")
        
        stop(m_err)
        
      }
      
    }
    
  }

  ########

  #### check ref_par ####
  
  if(!is.null(ref_par)){
  check_ref_par(ref_par = ref_par, parents = mppData$parents)
  }
  
  # Now split into GE and oneS

  if(GE){

    # 5.1 test format trait
    #################

    if(!is.character(trait)){

      stop("'trait' must be character")

    }

    if(length(trait) < 2){

      stop("'trait' must containt a least two trait values (one per environment)")

    }

    trait.names <- colnames(mppData$pheno)

    if (!all(trait %in% trait.names)){

      t_nms <- paste(trait.names, collapse = ', ')

      stop("'trait' must match the trait values from mppData: ", t_nms)

    }

    ##############

    # 5.2 test format VCOV
    ############

    if(!is.character(VCOV)){

      stop("VCOV must be character")

    }
    
    if(!fast){
      
      if (!(VCOV %in% c("ID", "CSRT", "CS_CSRT"))){
        
        stop("'VCOV' must be ", dQuote("ID"), ', ', dQuote("CSRT"), ' or ',
             dQuote("CS_CSRT"))
        
      }
      
    } else{
      
      if (!(VCOV %in% c("CS", "CSE", "CS_CSE", "UN", "UN_K"))){
        
        stop("'VCOV' must be one of: ", dQuote("CS"), ', ', dQuote("CSE"),
             dQuote("CS_CSE"),', ', dQuote("UN"), ' or ', dQuote("UN_K"))
        
      }
      
    }

    

    ###########

  } else { # one stage

    # 6.1 test format trait
    #################

    if(!is.character(trait)){

      stop("'trait' must be character")

    }

    trait.names <- colnames(plot_data)

    if (!all(trait %in% trait.names)){

      t_nms <- paste(trait.names, collapse = ', ')

      stop("'trait' must match the trait values from plot_data: ", t_nms)

    }

    ##############

    # 6.2 test format VCOV
    ############

    if(!is.character(VCOV)){

      stop("VCOV must be character")

    }

    if (!(VCOV %in% c("CSRT", "CS_CSRT", "CS_AR1xAR1", "CS_CSRT_AR1xAR1"))){

      stop("'VCOV' must be ", dQuote("CSRT"), ' or ',dQuote("CS_CSRT"))

    }

    ###########

    # 6.3 test format plot_data
    ###########

    if(!is.data.frame(plot_data)){

      stop("plot_data must be a data.frame")

    }

    pd_nms <- colnames(plot_data)

    if(!all(c('genotype', 'check', 'cross', 'env') %in% pd_nms)){

      stop("plot_data must contain the following columns: genotype, check, cross, and env")

    }

    # test if the plot data contain at least one value for each genotypes

    test_geno <- !all(mppData$geno.id %in% unique(plot_data$genotype))

    if(test_geno){

      prob_geno <- mppData$geno.id[!(mppData$geno.id %in% unique(plot_data$genotype))]

      message <- paste('The following genotypes:', paste(prob_geno, collapse = ', '),
                       'are present in the mppData object but do not have a value',
                       'in plot_data. Add a value for these genotypes in plot_data',
                       'or remove them from the mppData object.')

      stop(message)

    }


    ###########

    # 6.4 test format exp_des_form
    ###########

    if(!is.character(VCOV)){

      stop("exp_des_form must be character")

    }

    ###########

  }


}
