#############
# sign.star #
#############

# function that convert significance of results in stars

sign.star <- function(x){
  
  if(is.na(x)){
    
    sign <- NA
  
    } else {
      
      if(x>=0.1){
        sign <- ""
      }else if((0.1>x) & (0.05<=x)){
        sign <- "."
      }else if((0.05>x) & (0.01<=x)){
        sign <- "*"
      }else if((0.01>x) & (0.001<=x)){
        sign <- "**"
      }else if(0.001>x){
        sign <- "***"
      }
      
      return(sign)
      
    }
  
}