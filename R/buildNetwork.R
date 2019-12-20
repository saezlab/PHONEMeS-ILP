#'\code{buildNetwork}
#'
#'Network building from PKNlist object
#'

buildNetwork <- function(experiments=experiments, 
                         data.P=data.P, 
                         bg=bg, 
                         nK = "all", 
                         temp = temp){
  
  pkn <- build_Nw(data.On = data.P, targets.On = targets.P, bg = bg, nK = nK)
  
  pknList <- pkn
  
  interactions <- pknList@interactions
  species <- pknList@species
  
  interactions$SID <- NA
  
  for(i in 1:nrow(temp)){
    
    toAdd <- 
      as.data.frame(t(as.matrix(c(NA, NA, NA, temp[i, 1], 
                                  NA, NA, NA, temp[i, 3], 1))))
    colnames(toAdd) <- colnames(interactions)
    
    interactions <- rbind(interactions, toAdd)
    
  }
  
  tempSpecies <- unique(c(temp[, 1], temp[, 3]))
  
  idx1 <- grep(pattern = "_R1", x = tempSpecies)
  idx2 <- grep(pattern = "_S", x = tempSpecies)
  idx3 <- grep(pattern = "_T", x = tempSpecies)
  idx4 <- grep(pattern = "_Y", x = tempSpecies)
  sites <- tempSpecies[c(idx1, idx2, idx3, idx4)]
  
  for(i in 1:length(sites)){
    
    split <- strsplit(x = sites[i], split = "_", fixed = TRUE)[[1]]
    
    protein <- ""
    for(j in 1:(length(split)-1)){
      
      if(j==1){
        
        protein <- paste0(protein, split[j])
        
      }
      else{
        
        protein <- paste0(protein, "_", split[j])
        
      }
      
    }
    
    toAdd <- 
      as.data.frame(t(as.matrix(c(NA, NA, NA, sites[i], 
                                  NA, NA, NA, protein, 1))))
    colnames(toAdd) <- colnames(interactions)
    interactions <- rbind(interactions, toAdd)
    
  }
  
  dd <- interactions[, c(4, 8)]
  if(length(which(duplicated(dd))) > 0){
    
    interactions <- interactions[-which(duplicated(dd)), ]
    
  }
  
  interactions <- unique(interactions)
  
  interactions$SID[which(is.na(interactions$S.AC))] <- 
    paste0("i", 1:length(which(is.na(interactions$S.AC))))
  interactions$SID[which(!is.na(interactions$S.AC))] <- 
    paste0("e", 1:length(which(!is.na(interactions$S.AC))))
  
  pknList@interactions <- interactions
  pknList@species <- unique(c(pkn@species, temp[, 1], temp[, 3]))
  
  return(pknList)
  
}