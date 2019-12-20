#'\code{write_equality_constraints}
#'
#'Witig the equality constraints for each interaction variable
#'

write_equality_constraints <- function(dataMatrix = dataMatrix, 
                                       binaries = binaries, pknList = pknList){
  
  dM <- dataMatrix[[1]]
  equalityConstraints <- c()
  sif <- createSIF(pknList = pknList)
  
  for(ii in 1:nrow(dM)){
    
    for(jj in 1:nrow(sif)){
      
      aa <- binaries[[1]][which(
        binaries[[3]]==paste0("reaction ", sif[jj, 1], "=", sif[jj, 3]))]
      bb <- binaries[[1]][which(
        binaries[[3]]==paste0("interaction ", sif[jj, 1], "=", sif[jj, 3], 
                              " in experiment ", ii))]
      
      equalityConstraints <- c(equalityConstraints, 
                               paste0(aa, " - ", bb, " >= 0"))
      
    }
    
  }
  
  ll = list()
  for(ii in 1:nrow(sif)){
    
    ll[[length(ll)+1]] = ""
    for(jj in 1:nrow(dM)){
      
      ll[[ii]] = c(ll[[ii]],
                   binaries[[1]][which(binaries[[3]]==
                                         paste0("interaction ", sif[ii, 1], "=", 
                                                sif[ii, 3], " in experiment ", 
                                                jj))])
      
    }
    
    ll[[ii]] = ll[[ii]][-1]
    names(ll)[ii] = 
      binaries[[1]][which(binaries[[3]]==paste0("reaction ", sif[ii, 1], 
                                                "=", sif[ii, 3]))]
    
  }
  
  for(ii in 1:length(ll)){
    
    const = names(ll)[ii]
    
    for(jj in 1:length(ll[[ii]])){
      
      const = paste0(const, " - ", ll[[ii]][jj])
      
    }
    
    const = paste0(const, " <= 0")
    
    equalityConstraints = c(equalityConstraints, const)
    
  }
  
  return(equalityConstraints)
  
}