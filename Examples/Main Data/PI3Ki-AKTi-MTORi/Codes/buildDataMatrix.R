buildDataMatrix <- function(dataGMM = dataGMM, pknList = pknList, targets, experiments){
  
  tg <- unlist(targets)
  
  idxTG <- which(pknList@species==tg)
  idxDS <- which(!is.na(match(pknList@species, intersect(dataGMM@IDmap$S.cc, pknList@species))))
  idxDN <- setdiff(1:length(pknList@species), c(idxTG, idxDS))
  
  cNames <- c()
  for(i in 1:length(pknList@species)){
    
    if(i %in% idxTG){
      
      cNames <- c(cNames, paste0("TG:", pknList@species[i]))
      
    }
    
    if(i %in% idxDN){
      
      cNames <- c(cNames, paste0("DN:", pknList@species[i]))
      
    }
    
    if(i %in% idxDS){
      
      cNames <- c(cNames, paste0("DS:", pknList@species[i]))
      
    }
    
  }
  
  dataMatrix <- matrix(, nrow = length(experiments), ncol = length(pknList@species))
  colnames(dataMatrix) <- cNames
  dataMatrix[, c(idxTG, idxDN)] <- 0
  
  for(i in 1:length(experiments)){
    
    dM <- matrix(, nrow = length(experiments[[i]]), ncol = length(idxDS))
    colnames(dM) <- cNames[idxDS]
    rownames(dM) <- experiments[[i]]
    
    for(ii in 1:nrow(dM)){
      
      for(jj in 1:ncol(dM)){
        
        site <- strsplit(colnames(dM)[jj], split = ":")[[1]][2]
        siteID <- dataGMM@IDmap$dataID[which(dataGMM@IDmap$S.cc==site)][1]
        id <- which(names(dataGMM@res)==siteID)
        
        if(!is.na(dataGMM@res[[id]][rownames(dM)[ii], 4])){
          
          if(dataGMM@res[[id]][rownames(dM)[ii], 4]=="OK"){
            
            dM[ii, jj] <- as.numeric(dataGMM@res[[id]][rownames(dM)[ii], 1])
            
          }
          
        }
        
      }
      
    }
    
    for(j in 1:length(idxDS)){
      
      if(!all(is.na(dM[, j]))){
        
        dataMatrix[i, idxDS[j]] <- min(dM[, j], na.rm = TRUE)
        
      }
      
    }
    
  }
  
  dataMatrix[which(as.character(dataMatrix)=="-Inf")] <- -9999
  dataMatrix[which(as.character(dataMatrix)=="Inf")] <- 9999
  naIdx <- c()
  for(i in 1:length(idxDS)){
    
    if(all(is.na(dataMatrix[, idxDS[i]]))){
      
      naIdx <- c(naIdx, idxDS[i])
      
    }
    
  }
  
  if(length(naIdx) > 0){
    
    for(kk in 1:length(naIdx)){
      
      idxDN <- c(idxDN, naIdx[kk])
      idxDS <- idxDS[-which(idxDS==naIdx[kk])]
      
      cNames[naIdx[kk]] <- paste0("DN:", strsplit(cNames[naIdx[kk]], split = ":")[[1]][2])
      
    }
    
  }
  colnames(dataMatrix) <- cNames
  dataMatrix[is.na(dataMatrix)] <- 0
  
  res <- colSums(x = dataMatrix)
  
  res <- t(as.matrix(res))
  dataMatrix <- res
  
  res <- list(dataMatrix=dataMatrix, tgID=idxTG, dnID=idxDN, dsID=idxDS, species=pknList@species)
  
  return(res)
  
}
