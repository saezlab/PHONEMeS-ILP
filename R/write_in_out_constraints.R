#'\code{write_in_out_constraints}
#'
#'Writing the edge in-out constraints
#'

write_in_out_constraints <- function(binaries = binaries, 
                                     targets = targets, 
                                     dataMatrix = dataMatrix, 
                                     pknList = pknList){
  
  if((any(grepl(pattern = "orin", x = binaries[[3]]))) && 
     (any(grepl(pattern = "orout", x = binaries[[3]])))){
    
    constraints <- c()
    
    sif <- createSIF(pknList = pknList)
    
    #
    dnID <- c()
    for(ii in 1:nrow(dataMatrix[[1]])){dnID <- 
      c(dnID, dataMatrix$species[dataMatrix$dnID[[ii]]])}
    dsID <- c()
    for(ii in 1:nrow(dataMatrix[[1]])){dsID <- 
      c(dsID, dataMatrix$species[dataMatrix$dsID[[ii]]])}
    uSites <- c()
    for(ii in 1:length(dsID)){
      idx1 <- which(sif[, 1]==dsID[ii])
      idx2 <- which(sif[, 3]==dsID[ii])
      
      if((length(idx1)>0) && (length(idx2)>0)){uSites <- c(uSites, dsID[ii])}
      
    }
    dnID <- c(dnID, dsID)
    dnID <- unique(dnID)
    dnID <- setdiff(dnID, unique(unlist(targets)))
    
    for(ii in 1:length(dnID)){
      
      idx1 <- which(sif[, 3]==dnID[ii])
      idx2 <- which(sif[, 1]==dnID[ii])
      orinVar <- 
        binaries[[1]][which(binaries[[3]]==paste0("orin ", dnID[ii]))]
      oroutVar <- 
        binaries[[1]][which(binaries[[3]]==paste0("orout ", dnID[ii]))]
      
      if((length(idx1)>0) && (length(idx2)>0)){
        
        varIn <- c()
        for(jj in 1:length(idx1)){
          
          varIn <- 
            c(varIn, 
              binaries[[1]][which(binaries[[3]]==
                                    paste0("reaction ", sif[idx1[jj], 1], "=", 
                                           sif[idx1[jj], 3]))])
          
        }
        varOut <- c()
        for(jj in 1:length(idx2)){
          
          varOut <- 
            c(varOut, 
              binaries[[1]][which(binaries[[3]]==
                                    paste0("reaction ", sif[idx2[jj], 1], "=", 
                                           sif[idx2[jj], 3]))])
          
        }
        
        cc1 <- paste0(orinVar, " - ", oroutVar, " = 0")
        
        cc2 <- ""
        for(kk in 1:length(varIn)){
          if(kk==1){
            cc2 <- paste0(cc2, varIn[kk])
          } else {
            cc2 <- paste0(cc2, " + ", varIn[kk])
          }
        }
        cc2 <- paste0(cc2, " - ", orinVar, " >= 0")
        
        cc3 <- ""
        for(kk in 1:length(varOut)){
          if(kk==1){
            cc3 <- paste0(cc3, varOut[kk])
          } else {
            cc3 <- paste0(cc3, " + ", varOut[kk])
          }
        }
        cc3 <- paste0(cc3, " - ", oroutVar, " >= 0")
        
        cc4 <- c()
        for(kk in 1:length(varIn)){
          cc4 <- c(cc4, paste0(varIn[kk], " - ", orinVar, " <= 0"))
        }
        
        cc5 <- c()
        for(kk in 1:length(varOut)){
          cc5 <- c(cc5, paste0(varOut[kk], " - ", oroutVar, " <= 0"))
        }
        
        constraints <- c(constraints, c(cc1, cc2, cc3, cc4, cc5))
        
      }
      
    }
    
    return(constraints)
    
  } else {
    
    return(NULL)
    
  }
  
}
