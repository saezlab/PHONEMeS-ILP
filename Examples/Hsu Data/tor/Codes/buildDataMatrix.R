buildDataMatrix <- function(dataGMM = dataGMM, pknList = pknList, targets, experiments){
  
  tg <- unlist(targets)
  
  allSpecies <- pknList@species
  dataMatrix <- matrix(, nrow = length(unlist(experiments)), ncol = length(allSpecies))
  
  #colnames
  cNames <- c()
  kk <- which(allSpecies%in%tg)
  for(i in 1:length(kk)){
    cNames <- c(cNames, paste("TG:", allSpecies[kk[i]], sep = ""))
  }
  
  idxDN <- which(!(allSpecies%in%dataGMM@IDmap$S.cc))
  idxDN <- setdiff(idxDN, kk)
  for(i in 1:length(idxDN)){
    cNames <- c(cNames, paste("DN:", allSpecies[idxDN[i]], sep = ""))
  }
  
  dataMatrix[, 1:length(cNames)] <- as.numeric(0)
  
  idxDS <- which(allSpecies%in%dataGMM@IDmap$S.cc)
  ds <- length(cNames)
  for(i in 1:length(idxDS)){
    cNames <- c(cNames, paste("DS:", allSpecies[idxDS[i]], sep = ""))
  }
  
  colnames(dataMatrix) <- cNames
  
  #set values on the matrix
  for(i in 1:length(idxDS)){
    cpSite <- allSpecies[idxDS[i]]
    dataSite <- dataGMM@IDmap$dataID[which(dataGMM@IDmap$S.cc==cpSite)]
    for(j in 1:length(dataSite)){
      dataID <- which(names(dataGMM@res)==dataSite[j])
      for(dd in 1:length(dataID)){
        cnt <- 1
        for(k in 1:length(experiments)){
          for(l in 1:length(experiments[[k]])){
            if(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 4] == "OK"){
              if(is.na(dataMatrix[cnt, ds+i])){
                if(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 2] == "P"){
                  dataMatrix[cnt, ds+i] <- as.numeric(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 1])
                }
                else{
                  dataMatrix[cnt, ds+i] <- abs(as.numeric(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 1]))
                }
                #cnt <- cnt + 1
              }
              else{
                if(as.numeric(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 1]) < dataMatrix[cnt, ds+i]){
                  if(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 2] == "P"){
                    dataMatrix[cnt, ds+i] <- as.numeric(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 1])
                  }
                  else{
                    dataMatrix[cnt, ds+i] <- abs(as.numeric(dataGMM@res[[dataID[dd]]][which(rownames(dataGMM@res[[dataID[dd]]])==experiments[[k]][l]), 1]))
                  }
                }
              }
              cnt <- cnt + 1
            }
          }
        }
      }
    }
  }
  dataMatrix[is.na(dataMatrix)] <- as.numeric(0)
  dataMatrix[is.infinite(dataMatrix)] <- as.numeric(10000)
  
  #indeces
  tgID <- c()
  dnID <- c()
  dsID <- c()
  for(i in 1:ncol(dataMatrix)){
    
    if(strsplit(colnames(dataMatrix)[i], split = ":")[[1]][1]=="TG"){
      tgID <- c(tgID, i)
    }
    
    if(strsplit(colnames(dataMatrix)[i], split = ":")[[1]][1]=="DN"){
      dnID <- c(dnID, i)
    }
    
    if(strsplit(colnames(dataMatrix)[i], split = ":")[[1]][1]=="DS"){
      dsID <- c(dsID, i)
    }
    
  }
  
  #species
  species <- c()
  for(i in 1:ncol(dataMatrix)){
    
    species <- c(species, strsplit(colnames(dataMatrix)[i], split = ":")[[1]][2])
    
  }
  
  res <- list(dataMatrix=dataMatrix, tgID=tgID, dnID=dnID, dsID=dsID, species=species)
  
  return(res)
}