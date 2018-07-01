#################################
# check object format and class #
#################################

is_mppData <- function(x){
  
  inherits(x = x, what = "mppData")
  
}

# check trait
#############

check_trait <- function(trait, mppData){
  
  if(!(is.character(trait) || is.numeric(trait))){
    
    stop("'trait' must be numeric or character")
    
  }
  
  if(is.numeric(trait)){
    
    nb.trait <- dim(mppData$pheno)[2]
    
    if(!((0 < trait) && (trait <=nb.trait))){
      
      stop("trait must be between 1 and ", nb.trait,
           " the total number of traits in 'mppData'")
      
    }
    
  }
  
  if(is.character(trait)){
    
    trait.names <- colnames(mppData$pheno)
    
    if (!(trait %in% trait.names)){
      
      t_nms <- paste(trait.names, collapse = ', ')
      
      stop("'trait' must be one of: ", t_nms)
      
    }
    
  }
  
}

# peak the trait
################

sel_trait <- function(mppData, trait){
  
  if(is.numeric(trait)){
    
    t_val <- mppData$pheno[, trait]
    
  } else {
    
    trait.names <- colnames(mppData$pheno)
    t_val <- mppData$pheno[, which(trait.names == trait)]
    
  }
  
  return(t_val)
  
}

# check mppData
###############

# function to check that the mppData has a correct format
# check that if the user want to compute the ancestral model
# the parent clustering object is provided

check_mppData <- function(mppData, Q.eff = NULL){
  
  if(!is_mppData(mppData)) {
    
    stop("'mppData' must be of class ", dQuote("mppData"))
    
  }
  
  # check mppData is processed at least up to IBD
  
  if((mppData$status != 'IBD') && (mppData$status != 'complete')){
    
    stop("'mppData' is not complete. Use first all processing ",
         "functions in the specified order: QC.mppData, IBS.mppData, ",
         "IBD.mppData, and optionally parent_cluster.mppData for the ancestral ",
         "model")
    
  }
  
  if(!is.null(Q.eff)){
    
    if("anc" %in% Q.eff){
      
      if(mppData$status != 'complete'){
        
        stop("'mppData' do not contains the parent clustering information ",
             "necessary for the ancestral model.",
             "to add this information use parent_cluster.mppData")
        
      }
      
      
    }
    
  }
  
}
