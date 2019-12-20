#'\code{create_binary_variables_for_y_vector}
#'
#'Creating binary variables for the edges
#'

create_binary_variables_for_y_vector <- function(pknList = pknList){
  
  sif <- createSIF(pknList = pknList)
  
  y_vector <- paste0(sif[, 1], "=", sif[, 3])
  
  variables <- paste0("y_", 1:length(y_vector))
  
  binary_variables_list = list(1:(length(y_vector)),
                               variables, paste("reaction",y_vector))
  
  return(binary_variables_list)
  
}