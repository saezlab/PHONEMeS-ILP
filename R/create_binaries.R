#'\code{create_binaries}
#'
#'Creating the object containing the binary variables
#'

create_binaries <- function(binaries_x = binaries_x, 
                            binaries_z = binaries_z, 
                            binaries_in_out = binaries_in_out, 
                            binaries_y = binaries_y){
  
  numbers <- c()
  bins <- append(binaries_x[[1]], binaries_z[[1]]+length(binaries_x[[1]]))
  bins <- append(bins, binaries_in_out[[1]]+length(bins))
  bins <- append(bins, binaries_y[[1]]+length(bins))
  for(i in 1:length(bins)){
    numbers <- c(numbers, paste("xb", bins[i], sep = ""))
  }
  variables <- append(binaries_x[[2]], binaries_z[[2]])
  variables <- append(variables, binaries_in_out[[2]])
  variables <- append(variables, binaries_y[[2]])
  identifiers <- append(binaries_x[[3]], binaries_z[[3]])
  identifiers <- append(identifiers, binaries_in_out[[3]])
  identifiers <- append(identifiers, binaries_y[[3]])
  
  binaries = list(numbers, variables, identifiers)
  
  return(binaries)
  
}