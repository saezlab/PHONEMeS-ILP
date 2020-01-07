#'\code{write_constraints_2}
#'
#'Wting of constraints - 2
#'

write_constraints_2 <- function(dataMatrix = dataMatrix, 
                                binaries = binaries, 
                                pknList = pknList){
  
  constraints2 <- c()
  
  for(ii in 1:nrow(dataMatrix$dataMatrix)){
    
    tNames <- dataMatrix$species[dataMatrix$tgID[[ii]]]
    
    sif <- createSIF(pknList)
    
    
    tAdjacent <- list()
    for(i in 1:length(tNames)){
      
      temp <- c()
      tAdjacent[[length(tAdjacent)+1]] <- c(temp, which(sif[, 1]==tNames[i]))
      
    }
    
    for(i in 1:length(tAdjacent)){
      
      temp <- ""
      
      if(length(tAdjacent[[i]])>0){
        
        for(j in 1:length(tAdjacent[[i]])){
          
          if(j==1){
            
            temp <- 
              binaries[[1]][which(binaries[[3]]==
                                    paste0("interaction ", 
                                           sif[tAdjacent[[i]][j], 1], "=", 
                                           sif[tAdjacent[[i]][j], 3], 
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
        
        constraints2 <- c(constraints2, paste(temp, " >= 1", sep = ""))
        
      }
      
    }
    
    for(i in 1:length(tNames)){
      
      constraints2 <- 
        c(constraints2, 
          paste0(binaries[[1]][which(binaries[[3]]==
                                       paste0("species ", tNames[i], 
                                              " in experiment ", ii))], " = 1"))
      
    }
    
  }
  
  return(constraints2)
  
}