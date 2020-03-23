#'\code{write_objective_function}
#'
#'Creating the object containing the binary variables
#'

write_objective_function <- function(dataMatrix = dataMatrix, 
                                     binaries = binaries, sizePen = 0.0001){
  
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
    
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "interaction"){
      
      # r1 <- strsplit(
      #   strsplit(binaries[[3]][i], split = " ")[[1]][2], split = "=")[[1]][1]
      # r2 <- strsplit(
      #   strsplit(binaries[[3]][i], split = " ")[[1]][2], split = "=")[[1]][2]
      # 
      # if(grepl(pattern = "_R1", x = r1) || grepl(pattern = "_R1", x = r2)){
      #   
      #   cc <- c(grepl(pattern = "_R1", x = r1), grepl(pattern = "_R1", x = r2))
      #   if(cc[1]==FALSE && cc[2]==TRUE){
      #     
      #     objectiveFunction <- 
      #       paste(objectiveFunction, " - 0.001 ", binaries[[1]][i], sep = "")
      #     
      #   }
      #   else{
      #     
      #     objectiveFunction <- 
      #       paste(objectiveFunction, " - 0.1 ", binaries[[1]][i], sep = "")
      #     
      #   }
      #   
      # }
      # else{
      #   
      #   objectiveFunction <- paste(objectiveFunction, " - ", 
      #                              sizePen, " ", binaries[[1]][i], sep = "")
      #   
      # }
      
      objectiveFunction <- paste(objectiveFunction, " + ", 
                                 sizePen, " ", binaries[[1]][i], sep = "")
      
    }
    
  }
  
  return(objectiveFunction)
  
}
