###################
# reference.count #
###################

# counter to follow the state of permutation test

reference.count <- function(N, l=10){
  
  # decompose the value
  
  if(N >= 100) {
    
    limit <- seq(0,N,(N/l))
    limit <- limit[-1]
    percentage <- seq(0,100,by=l)
    percentage <- percentage[-1]
    
  }
  
  return(data.frame(limit,percentage))
  
}

reference.count(N = 2000)