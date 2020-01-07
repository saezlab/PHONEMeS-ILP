#'\code{create_orin_orout_variables}
#'
#'Creating the or-in-or-out (gamma) variables
#'

create_orin_orout_variables <- function(dataMatrix = dataMatrix, 
                                        targets = targets, pknList = pknList){
  
  dM = dataMatrix[[1]]
  
  if(nrow(dM) > 1){
    
    dnID <- c()
    for(ii in 1:nrow(dataMatrix[[1]])){dnID <- 
      c(dnID, dataMatrix$species[dataMatrix$dnID[[ii]]])}
    dsID <- c()
    for(ii in 1:nrow(dataMatrix[[1]])){dsID <- 
      c(dsID, dataMatrix$species[dataMatrix$dsID[[ii]]])}
    uSites <- c()
    sif <- createSIF(pknList = pknList)
    for(ii in 1:length(dsID)){
      idx1 <- which(sif[, 1]==dsID[ii])
      idx2 <- which(sif[, 3]==dsID[ii])
      
      if((length(idx1)>0) && (length(idx2)>0)){uSites <- c(uSites, dsID[ii])}
      
    }
    dnID <- c(dnID, dsID)
    dnID <- unique(dnID)
    dnID <- setdiff(dnID, unique(unlist(targets)))
    
    variables <- paste0("orin_", 1:length(dnID))
    variables <- c(variables, paste0("orout_", 1:length(dnID)))
    
    identifiers <- paste0("orin ", dnID)
    identifiers <- c(identifiers, paste0("orout ", dnID))
    
    binary_variables_list = list(1:length(variables),variables, identifiers)
    
    return(binary_variables_list)
    
  } else {
    
    return(NULL)
    
  }
  
}