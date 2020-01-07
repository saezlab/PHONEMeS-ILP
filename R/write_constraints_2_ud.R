#'\code{write_constraints_2_ud}
#'
#'Writing the constraints - 2 for the upside-down case
#'

write_constraints_2_ud <- function(dataMatrix = dataMatrix, 
                                   binaries = binaries, 
                                   pknList = pknList){
  
  constraints2 <- c()
  
  for(ii in 1:nrow(dataMatrix$dataMatrix)){
    
    tNames <- dataMatrix$species[dataMatrix$tgID[[ii]]]
    
    sif <- createSIF(pknList)
    
    tAdjacent <- list()
    for(i in 1:length(tNames)){
      
      temp <- c()
      tAdjacent[[length(tAdjacent)+1]] <- c(temp, which(sif[, 3]==tNames[i]))
      
    }
    
    for(i in 1:length(tAdjacent)){
      
      temp <- ""
      
      if(length(tAdjacent[[i]])>0){
        
        for(j in 1:length(tAdjacent[[i]])){
          
          if(j==1){
            
            temp <- 
              binaries[[1]][which(binaries[[3]]==
                                    paste0("interaction ", 
                                           sif[tAdjacent[[i]][j], 1], 
                                           "=", sif[tAdjacent[[i]][j], 3], 
                                           " in experiment ", ii))]
            
          }
          else{
            
            temp <- 
              paste0(temp, " + ", 
                     binaries[[1]][which(binaries[[3]]==
                                           paste0("interaction ", 
                                                  sif[tAdjacent[[i]][j], 1], 
                                                  "=", 
                                                  sif[tAdjacent[[i]][j], 3], 
                                                  " in experiment ", ii))])
            
          }
          
        }
        
        constraints2 <- 
          c(constraints2, 
            paste(temp, " - ", 
                  binaries[[1]][which(binaries[[3]]==
                                        paste0("species ", tNames[i], 
                                               " in experiment ", ii))], 
                  " >= 0", sep = ""))
        
      }
      
    }
    
  }
  
  return(constraints2)
  
}