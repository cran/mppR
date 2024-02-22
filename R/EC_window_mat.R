#################
# EC_window_mat #
#################

# Function to determine all possible windows or 10, 20, ... x days given
# the crop duration (longest window = crop_duration - 5 days), and
# starting dates.

EC_window_mat <- function(crop_duration, min_win = 20){
  
  windows <- seq(0, crop_duration, min_win)
  windows <- windows[-which(windows == 0)]
  n_win <- length(windows)
  start_day <- 1:(crop_duration - min_win)
  
  # for each cycle break we can calculate a certain number of windows.
  end_day_mat <- matrix(NA, nrow = n_win, ncol = length(start_day))
  
  for(i in 1:n_win){
    
    test <- TRUE
    st_d <- 1
    
    while(test){
      
      end_day <- windows[i] + start_day[st_d] - 1
      
      if(!(end_day > crop_duration)){
        
        end_day_mat[i, st_d] <- end_day
        st_d <- st_d + 1
        
        if(st_d > length(start_day)) {test <- FALSE}
        
      } else {test <- FALSE}
      
    }
    
  }
  
  rownames(end_day_mat) <- windows
  
  return(end_day_mat)
  
}