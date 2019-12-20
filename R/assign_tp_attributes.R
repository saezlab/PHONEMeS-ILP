#'\code{assign_tp_attributes}
#'
#'Assigning the time-point attributes for each edge present in the final 
#'solution
#'

assign_tp_attributes <- function(sifList = sifList){
  
  returnSIF <- matrix(data = c("", 0, "", ""), nrow = 1, ncol = 4)
  colnames(returnSIF) <- c("Source", "f50", "tp", "Target")
  
  for(ii in 1:length(sifList[[1]])){
    
    for(jj in 1:length(sifList)){
      
      currSIF <- sifList[[jj]][[ii]]
      
      for(kk in 1:nrow(currSIF)){
        
        ss <- currSIF[kk, 1]
        tt <- currSIF[kk, 3]
        
        idx1 <- which(returnSIF[, 1]==ss)
        idx2 <- which(returnSIF[, 4]==tt)
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx)<=0){
          
          returnSIF <- rbind(returnSIF, 
                             c(ss, "1", paste0("tp_", as.character(ii)), tt))
          
        } else {
          
          returnSIF[idx, 2] <- as.character(as.numeric(returnSIF[idx, 2])+1)
          
        }
        
      }
      
    }
    
  }
  
  returnSIF <- returnSIF[-1, ]
  
  for(ii in 1:length(sifList[[1]])){
    
    returnSIF[which(returnSIF[, 3]==paste0("tp_", ii)), 2] <- 
      as.character(
        as.numeric(
          returnSIF[which(returnSIF[, 3]==
                            paste0("tp_", ii)), 2])/(length(sifList[[1]])-ii+1))
    
  }
  
  return(returnSIF)
  
}