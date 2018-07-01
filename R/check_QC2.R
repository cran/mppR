############
# check_QC2 #
############

# check the data format for the quality control procedure.



check_QC2 <- function(mppData, n.lim, MAF.pop.lim, mk.miss, gen.miss,
                      MAF.cr.lim, MAF.cr.lim2, n.cores){
  
  # test the format of the mppData
  
  if(!is_mppData(mppData)){
    
    stop("'mppData' must be of class ", dQuote("mppData"))
    
  }
  
  # test if correct step in the mppData processing
  
  if(mppData$status != 'init'){
    
    stop("you have to process 'mppData' in a strict order: ",
         "create.mppData, QC.mppData, IBS.mppData, IBD.mppData, ",
         "parent_cluster.mppData. You can only use QC.mppData ",
         "after create.mppData")
    
  }
  
  # n.lim
  #######
  
  if(all(n.lim > table(mppData$cross.ind))){
    
    stop("the number of individual per cross is below 'n.lim' in all crosses")
    
  }
  
  # MAF.pop.lim
  #############
  
  # test that the values of MAF.pop.lim, mk.miss and gen.miss are
  # between 0 and 1.
  
  t1 <- ((MAF.pop.lim > 0) & (MAF.pop.lim < 1))
  
  if(!t1){
    
    stop("'MAF.pop.lim' can only take values between 0 and 1")
    
  }
  
  # mk.miss and gen.miss
  ######################
  
  t2 <- ((mk.miss >= 0) & (mk.miss <= 1))
  
  if(!t2){
    
    stop("'mk.miss' can only take values between 0 and 1")
    
  }
  
  t3 <- ((gen.miss >= 0) & (gen.miss <= 1))
  
  if(!t3){
    
    stop("'gen.miss' can only take values between 0 and 1")
    
  }
  
  # MAF.cr.lim
  ############
  
  # test if there are as many within cross MAF as cross in the MAF.cr.lim argument
  # test if the values are between 0 and 1.
  
  if(!is.null(MAF.cr.lim)){
    
    if(length(MAF.cr.lim) != mppData$n.cr){
      
      stop("'MAF.cr.lim' must contain one value per cross")
      
    }
    
    test <- any((MAF.cr.lim > 0) & (MAF.cr.lim < 1))
    
    if(test){
      
      stop("'MAF.cr.lim' can only contain values between 0 and 1")
      
    }
    
  }
  
  # MAF.cr.lim2
  #############
  
  if(!is.null(MAF.cr.lim2)){
    
    if(length(MAF.cr.lim2) !=1 ){
      
      stop("'MAF.cr.lim2' can only contain a single value between 0 and 1")
      
    }
    
    if(!((MAF.cr.lim2 > 0) & (MAF.cr.lim2 < 1))) {
      
      stop("'MAF.cr.lim2' can only contain a value between 0 and 1")
      
    }
    
  }
  
  # n.cores
  #########
  
  if(!is.numeric(n.cores)){
    
    stop("'n.cores must be numeric")
    
  }

}