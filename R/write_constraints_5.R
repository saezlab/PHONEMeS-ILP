#'\code{write_constraints_5}
#'
#'Wting of constraints - 5
#'

write_constraints_5 <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  constraints4 <- c()
  
  sif <- createSIF(pknList)
  
  for(ii in 1:nrow(dataMatrix$dataMatrix)){
    
    if(nrow(dataMatrix$dataMatrix)>1){
      
      nNames <- dataMatrix$species[dataMatrix$dsID[[ii]]]
      
    } else {
      
      nNames <- dataMatrix$species[dataMatrix$dsID[[1]]]
      
    }
    
    nIncident <- list()
    for(i in 1:length(nNames)){
      
      temp <- c()
      nIncident[[length(nIncident)+1]] <- c(temp, which(sif[, 3]==nNames[i]))
      
    }
    
    for(i in 1:length(nIncident)){
      if(length(nIncident[[i]]) > 0){
        temp <- ""
        for(j in 1:length(nIncident[[i]])){
          
          ss <- sif[nIncident[[i]][j], 1]
          tt <- sif[nIncident[[i]][j], 3]
          
          if(j==1){
            temp <- binaries[[1]][which(binaries[[3]]==paste0("interaction ", ss, "=", tt, " in experiment ", ii))]
          }
          else{
            temp <- paste0(temp, " + ", binaries[[1]][which(binaries[[3]]==paste0("interaction ", ss, "=", tt, " in experiment ", ii))])
          }
          
        }
        
        bb <- which(binaries[[3]]==paste0("species ", nNames[i], " in experiment ", ii))
        constraints4 <- c(constraints4, paste(temp, " - ", binaries[[1]][bb], " >= 0"))
      }
      else{
        
        bb <- which(binaries[[3]]==paste0("species ", nNames[i], " in experiment ", ii))
        constraints4 <- c(constraints4, paste0(binaries[[1]][bb], " = 0"))
        
      }
      
    }
    
  }
  
  if(length(grep(pattern = "NA", x = constraints4))==0){
    
    return(constraints4)
    
  }
  else{
    
    return(constraints4[-grep(pattern = "NA", x = constraints4)])
    
  }
  
}
