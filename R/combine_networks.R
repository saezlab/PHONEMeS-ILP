#'\code{combine_networks}
#'
#'Cobining the time-point networks
#'

combine_networks <- function(resList = resList){
  
  returnSIF <- resList[[1]]
  
  if(length(resList)==1){
    
    return(returnSIF)
    
  } else {
    
    for(ii in 2:length(resList)){
      
      currSIF <- resList[[ii]]
      
      for(jj in 1:nrow(currSIF)){
        
        idx1 <- which(returnSIF[, 1]==currSIF[jj, 1])
        idx2 <- which(returnSIF[, 3]==currSIF[jj, 3])
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx)>0){
          
          returnSIF[idx, 2] <- returnSIF[idx, 2] + 1
          
        } else {
          
          returnSIF <- rbind(returnSIF, currSIF[jj, ])
          
        }
        
      }
      
    }
    
  }
  
  return(returnSIF)
  
}