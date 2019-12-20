#'\code{create_binary_variables_for_z_vector}
#'
#'Creating binary variables for interactions on each experimental conditions
#'

create_binary_variables_for_z_vector <- function(pknList = pknList, 
                                                 dataMatrix = dataMatrix){
  
  sif <- createSIF(pknList = pknList)
  
  z_vector <- paste0(sif[, 1], "=", sif[, 3])
  
  variables <- c()
  for(ii in 1:nrow(dataMatrix$dataMatrix)){
    for(jj in 1:nrow(sif)){
      variables <- c(variables, paste0("z_",jj,"^",ii))
    }
  }
  
  identifiers <- c()
  for(ii in 1:nrow(dataMatrix$dataMatrix)){
    identifiers <- c(identifiers, paste0("interaction ", 
                                         z_vector, " in experiment ", ii))
  }
  
  binary_variables_list = list(1:length(variables),variables, identifiers)
  
  return(binary_variables_list)
}