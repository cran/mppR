#############
# check_IBS #
#############

# Function to check the argument provided to the function IBS.mppData which
# convert the genotype data into 0, 1, 2 format and do an optional
# genotype imputation.

# check_IBS <- function(mppData, impute, impute.type, map_bp, replace.value){

check_IBS <- function(mppData){
  
  # check mppData format
  
  if(!is_mppData(mppData)){
    
    stop("'mppData' must be of class ", dQuote("mppData"))
    
  }
  
  # test if correct step in the mppData processing
  
  if(mppData$status != 'QC'){
    
    stop("you have to process 'mppData' in a strict order: ",
         "create.mppData, QC.mppData, IBS.mppData, IBD.mppData, ",
         "parent_cluster.mppData. You can only use IBD.mppData ",
         "after create.mppData, QC.mppData, and IBS.mppData")
    
  }
  
  # if(impute){
  #   
  #   if(!(impute.type %in% c("random","family","beagle","beagleAfterFamily",
  #                           "beagleNoRand", "beagleAfterFamilyNoRand","fix"))){
  #     
  #     stop("'impute.type' must be ", dQuote("random"), ', ', dQuote("family"), ', ',
  #          dQuote("beagle"), ', ', dQuote("beagleAfterFamily"), ', ',
  #          dQuote("beagleNoRand"), dQuote("beagleAfterFamilyNoRand"),
  #          'or ', dQuote("fix"))
  #     
  #   }
  #   
  #   # check if impute is beagle the right format of the map_bp
  #   
  #   if (impute.type %in% c("beagle","beagleAfterFamily",
  #                          "beagleNoRand", "beagleAfterFamilyNoRand")){
  #     
  #     if(!("package:synbreed" %in% search())){
  #       
  #       stop("to make imputation using Beagle please load package synbreed")
  #       
  #     }
  #     
  #     if(is.null(map_bp)){
  #       
  #       stop("'map_pb' is not provided")
  #       
  #     }
  #     
  #     if(!identical(map_bp[, 1], colnames(mppData$geno.off))){
  #       
  #       stop("the marker indicators of 'map_bp' and 'geno.off' are not identical")
  #       
  #     }
  #     
  #     if(!(is.character(map_bp[, 2]) || is.numeric(map_bp[, 2]))){
  #       
  #       stop("the chromosome indicator of 'map_bp' must be character or numeric")
  #       
  #     }
  #     
  #   }
  #   
  #   
  #   # if impute.type is fix. Chekc that a value was provided for replace.value
  #   
  #   if((impute.type == "fix") & (is.null(replace.value))){
  #     
  #     stop("'replace.value' is not provided")
  #     
  #   }
  #   
  #   if((impute.type == "fix")){
  #     
  #     if(!is.numeric(replace.value)){
  #       
  #       stop("'replace.value' must be numeric")
  #       
  #     }
  #     
  #     if (!(replace.value %in% c(0, 1, 2))){
  #       
  #       stop("'replace.value' must be 0, 1 or 2")
  #       
  #     }
  #     
  #   }
  #   
  # }

  
}