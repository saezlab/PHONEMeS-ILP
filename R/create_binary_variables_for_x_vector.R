#'\code{create_binary_variables_for_x_vector}
#'
#'Creating SIF file from pknList object
#'

create_binary_variables_for_x_vector <- function(dataMatrix = dataMatrix){
  
  allSpecies = c()
  for(i in 1:ncol(dataMatrix[[1]])){
    allSpecies <- c(allSpecies, 
                    strsplit(colnames(dataMatrix[[1]])[i], ":")[[1]][2])
  }
  
  variables = c()
  for(i in 1:nrow(dataMatrix[[1]])){
    for(j in 1:length(allSpecies)){
      variables <- c(variables, paste("x_", j, "^", i, sep = ""))
    }
  }
  
  identifiers = c()
  for(i in 1:nrow(dataMatrix[[1]])){
    for(j in 1:length(allSpecies)){
      identifiers <- c(identifiers, 
                       paste("species", allSpecies[j], "in experiment", i))
    }
  }
  
  binary_variables_list = list(1:length(variables),variables, identifiers)
  
  return(binary_variables_list)
  
}