#'\code{readOutSIFAll}
#'
#'Reading out SIF file from all the retrieved CPLEX solutions
#'

readOutSIFAll<- function(cplexSolutionFileName, binaries = binaries){
  
  reacIDX <- c()
  for(i in 1:length(binaries[[3]])){
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
      reacIDX <- c(reacIDX, i)
    }
  }
  reacVar <- binaries[[1]][reacIDX]
  
  cplexSolutionData <- xmlParse(cplexSolutionFileName)
  cplexSolution <- xmlToList(cplexSolutionData)
  
  cplexSolutionEdgesAll <- list()
  sifAll <- list()
  
  for(ii in 1:(length(cplexSolution)-1)){
    
    sif <- matrix(data = NA, nrow = 1, ncol = 3)
    
    currSolution <- cplexSolution[[ii]][[4]]
    
    for(jj in 1:length(currSolution)){
      
      if((currSolution[[jj]][1]%in%reacVar) && 
         round(as.numeric(currSolution[[jj]][3])==1)){
        
        reaction <- 
          strsplit(binaries[[3]][which(binaries[[1]]==currSolution[[jj]][1])], 
                   split = " ")[[1]][2]
        
        sif2bind <- 
          as.matrix(c(strsplit(x = reaction, split = "=")[[1]][1], 
                      "1", strsplit(x = reaction, split = "=")[[1]][2]))
        sif2bind <- t(sif2bind)
        sif <- rbind(sif, sif2bind)
        
      }
      
    }
    
    sifAll[[length(sifAll)+1]] <- unique(sif[-1, ])
    
  }
  
  for(ii in 1:length(sifAll)){
    
    if(ii==1){
      
      sif <- sifAll[[1]]
      
    } else {
      
      # sif <- unique(rbind(sif, sifAll[[ii]]))
      if(nrow(sifAll[[ii]]) > 0){
        
        for(jj in 1:nrow(sifAll[[ii]])){
          
          idx1 <- which(sif[, 1]==sifAll[[ii]][jj, 1])
          idx2 <- which(sif[, 3]==sifAll[[ii]][jj, 3])
          
          idx <- intersect(x = idx1, y = idx2)
          
          if(length(idx) > 0){
            
            sif[idx, 2] <- as.character(as.numeric(sif[idx, 2])+1)
            
          } else {
            
            sif <- rbind(sif, sifAll[[ii]][jj, ])
            
          }
          
        }
        
      }
      
    }
    
  }
  
  return(sif)
  
}