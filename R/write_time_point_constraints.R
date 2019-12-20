#'\code{write_time_point_constraints}
#'
#'Writing the time-point constraints
#'

write_time_point_constraints <- function(binaries = binaries, 
                                         tempSIF = tempSIF, 
                                         resultsSIF1 = resultsSIF1){
  
  idx <- c()
  for(i in 1:nrow(tempSIF)){
    
    kk <- 
      intersect(which(resultsSIF1[, 1]==tempSIF[i, 1]), 
                which(resultsSIF1[, 3]==tempSIF[i, 3]))
    
    if(length(kk) > 0){
      
      idx <- c(idx, kk)
      
    }
    
  }
  
  if(length(idx) > 0){
    
    constraints <- c()
    for(i in 1:length(idx)){
      
      cc <- 
        paste0(binaries[[1]][which(binaries[[3]]==
                                     paste0("reaction ", resultsSIF1[idx[i], 1], 
                                            "=", resultsSIF1[idx[i], 3]))])
      
      constraints <- c(constraints, paste0(cc, " = 1"))
      
    }
    
    return(constraints)
    
  }
  else{
    
    return(NULL)
    
  }
  
}