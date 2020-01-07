#'\code{removeRedundantNodes}
#'
#'Removing the redundant nodes from the directed protein-protein interactions
#'

removeRedundantNodes <- function(resultsSIF1 = resultsSIF1){
  
  idxToRem <- c()
  
  temp <- resultsSIF1
  
  for(i in 1:nrow(resultsSIF1)){
    
    temp[i, 1] <- gsub(pattern = "_R1", replacement = "", x = resultsSIF1[i, 1])
    temp[i, 3] <- gsub(pattern = "_R1", replacement = "", x = resultsSIF1[i, 3])
    
  }
  
  for(i in 1:nrow(resultsSIF1)){
    
    if(temp[i, 1]==temp[i, 3]){
      
      idxToRem <- c(idxToRem, i)
      
    }
    
  }
  
  return(temp[-idxToRem, ])
  
}