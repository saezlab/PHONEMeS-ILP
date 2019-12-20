#'\code{readOutSIF}
#'
#'Reading out SIF file from the CPLEX solution
#'

readOutSIF<- function(cplexSolutionFileName, binaries = binaries){
  
  reacIDX <- c()
  for(i in 1:length(binaries[[3]])){
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
      reacIDX <- c(reacIDX, i)
    }
  }
  reacVar <- binaries[[1]][reacIDX]
  
  cplexSolutionData <- xmlParse(cplexSolutionFileName)
  cplexSolution <- xmlToList(cplexSolutionData)
  
  cplexSolutionEdges <- list()
  for(i in 1:length(cplexSolution$CPLEXSolution$variables)){
    
    if(cplexSolution$CPLEXSolution$variables[i]$variable[1] %in% reacVar){
      
      cplexSolutionEdges[[length(cplexSolutionEdges)+1]] <- 
        cplexSolution$CPLEXSolution$variables[i]
      
    }
    
  }
  
  reac <- c()
  for(i in 1:length(cplexSolutionEdges)){
    
    if(round(as.numeric(cplexSolutionEdges[[i]]$variable[3])) == 1){
      
      reacTemp <- 
        binaries[[3]][which(binaries[[1]]==cplexSolutionEdges[[i]]$variable[1])]
      reac <- c(reac, reacTemp)
      
    }
    
  }
  
  sif <- matrix(, nrow = length(reac), ncol = 3)
  sif[, 2] <- 1
  for(i in 1:length(reac)){
    
    sif[i, 1] <- 
      strsplit(strsplit(reac[i], split = " ")[[1]][2], split = "=")[[1]][1]
    sif[i, 3] <- 
      strsplit(strsplit(reac[i], split = " ")[[1]][2], split = "=")[[1]][2]
    
  }
  
  return(sif)
  
}