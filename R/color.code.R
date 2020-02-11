##############
# color.code #
##############

# color indicator function

# all p-values from a nested model are transformed in a value between
# -6 and 6 which will give a indication of colour (blue or red) and the
# intensity of the colour

# color.code <- function(x){
#   
#   # red colours
#   
#   if(x<0){
#     
#     x <- abs(x)
#     
#     if(x>=0.05){
#       col <- 0
#     }else if((0.05>x) & (1e-2<=x)){
#       col <- -1
#     }else if((1e-2>x) & (1e-3<=x)){
#       col <- -2
#     }else if((1e-3>x) & (1e-4<=x)){
#       col <- -3
#     }else if((1e-4>x) & (1e-5<=x)){
#       col <- -4
#     }else if(1e-5>x){
#       col <- -5
#     }
#     
#   }else{
#     
#     # blue colors
#     
#     if(x>=0.05){
#       col <- 0
#     }else if((0.05>x) & (1e-2<=x)){
#       col <- 1
#     }else if((1e-2>x) & (1e-3<=x)){
#       col <- 2
#     }else if((1e-3>x) & (1e-4<=x)){
#       col <- 3
#     }else if((1e-4>x) & (1e-5<=x)){
#       col <- 4
#     }else if(1e-5>x){
#       col <- 5
#     }
#     
#   }
#   
#   return(col)
#   
# }

color.code <- function(x){
  
  # red colours
  
  if(x<0){
    
    x <- abs(x)
    
    if(x>=0.05){
      col <- 0
    }else if((0.05>x) & (1e-6<=x)){
      col <- -(-log10(x))
    }else if(1e-6>x){
      col <- -6
    }
    
  }else{
    
    # blue colors
    
    if(x>=0.05){
      col <- 0
    }else if((0.05>x) & (1e-6<=x)){
      col <- (-log10(x))
    }else if(1e-6>x){
      col <- 6
    }
    
  }
  
  return(col)
  
}