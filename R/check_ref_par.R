#################
# check_ref_par #
#################

# function to check the

check_ref_par <- function(ref_par, parents){
  
  if(!is.null(ref_par)){
    
    # test that there is only one reference parent and one connected part.
    
    if(length(ref_par) !=1){
      
      stop("You can only specify one 'ref_par'")
      
    }
    
    # test that reference parent is present in the list of parents.
    
    if(!(ref_par %in% parents)){
      
      par_list <- paste(parents, collapse = ", ")
      
      stop("'ref_par' must be one of: ", par_list)
      
    }
    
  }
  
}