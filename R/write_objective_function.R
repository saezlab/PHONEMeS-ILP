#'\code{write_objective_function}
#'
#'Creating the object containing the binary variables
#'

write_objective_function <- function(dataMatrix = dataMatrix, 
                                     binaries = binaries, sizePen = 0.00001){
  
  options(scipen=999)
  
  dM <- dataMatrix[[1]]
  kk <- 1
  objectiveFunction <- "obj:\t"
  for(i in 1:nrow(dM)){
    for(j in 1:ncol(dM)){
      if(as.numeric(dM[i, j]) <= 0){
        if(as.numeric(dM[i, j]) == 0){
          objectiveFunction <- paste(objectiveFunction, "", sep = "")
          kk <- kk+1
        }
        else{
          objectiveFunction <- 
            paste(objectiveFunction, " ", 
                  as.numeric(dM[i, j]), " ", "xb", kk, sep = "")
          kk <- kk+1
        }
      }
      else{
        objectiveFunction <- 
          paste(objectiveFunction, " + ", 
                as.numeric(dM[i, j]), " ", "xb", kk, sep = "")
        kk <- kk+1
      }
    }
  }
  
  for(i in 1:length(binaries[[3]])){
    
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
      
      objectiveFunction <- paste(objectiveFunction, " + ", 
                                 sizePen, " ", binaries[[1]][i], sep = "")
      
    }
    
  }
  
  return(objectiveFunction)
  
}