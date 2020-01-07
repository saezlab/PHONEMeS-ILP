#'\code{removeRedundantEdges}
#'
#'Removing the edges between the redundant nodes from the directed protein-
#'protein interactions
#'

removeRedundantEdges <- function(resultsSIF1 = resultsSIF1){
  
  kinases <- c()
  for(i in 1:nrow(resultsSIF1)){
    
    if(length(strsplit(resultsSIF1[i, 3], split = "_")[[1]])==2){
      
      kinases <- c(kinases, resultsSIF1[i, 3])
      
    }
    
  }
  kinases <- unique(kinases)
  
  toRem <- c()
  for(i in 1:length(kinases)){
    
    idx <- which(resultsSIF1[, 3]==kinases[i])
    
    if(length(idx) > 1){
      
      ctrl <- 0
      for(j in 1:length(idx)){
        
        if(strsplit(resultsSIF1[idx[j], 1], split = "_")[[1]][3] != "R1"){
          
          ctrl <- ctrl +1
          
        }
        
      }
      
      if(ctrl > 0){
        
        for(j in 1:length(idx)){
          
          if(strsplit(resultsSIF1[idx[j], 1], split = "_")[[1]][3] == "R1"){
            
            toRem <- c(toRem, idx[j])
            
          }
          
        }
        
      }
      
    }
    
  }
  
  if(length(toRem) > 0){
    
    idxRem <- c()
    
    speciesRem <- resultsSIF1[toRem, 1]
    
    for(i in 1:nrow(resultsSIF1)){
      
      if(resultsSIF1[i, 1]%in%speciesRem || resultsSIF1[i, 3]%in%speciesRem){
        
        idxRem <- c(idxRem, i)
        
      }
      
    }
    
    resultsSIF1 <- resultsSIF1[-idxRem, ]
    
    return(resultsSIF1)
    
  }
  else{
    
    return(resultsSIF1)
    
  }
  
}