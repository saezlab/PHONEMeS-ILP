#
#  This file is part of the CNO software
#
#  Copyright (c) 2018 - RWTH Aachen - JRC COMBINE
#
#  File author(s): E. Gjerga (enio.gjerga@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: https://saezlab.github.io/PHONEMeS/
#
##############################################################################
# $Id$

# Building data-matrix function which contains measurements & scores

buildDataMatrix <- function(dataGMM = dataGMM, pknList = pknList, targets, experiments){
  
  TG <- unlist(unique(targets))
  sif = createSIF(pknList = pknList)
  if(length(setdiff(x = TG, unique(c(sif[, 1], sif[, 3]))))>0){TG <- TG[-which(TG %in% setdiff(x = TG, unique(c(sif[, 1], sif[, 3]))))]}
  if(length(which(TG==""))>0){TG <- TG[-which(TG=="")]}
  sites <- intersect(pknList@species, dataGMM@IDmap$S.cc)
  idxTG <- list()
  idxDN <- list()
  idxDS <- list()
  
  dataMatrix <- matrix(data = 0, nrow = length(experiments), ncol = length(pknList@species))
  cNames <- unique(c(TG, setdiff(pknList@species, c(TG, sites)), sites))
  if(length(setdiff(x = cNames, y = pknList@species))>0){cNames <- cNames[-which(cNames==setdiff(x = cNames, y = pknList@species))]}
  idx <- which(cNames%in%sites)
  for(ii in 1:length(experiments)){
    idxDS[[length(idxDS)+1]] <- idx
  }
  cNames[setdiff(1:length(cNames), idx)] <- paste0("DN:", cNames[setdiff(1:length(cNames), idx)])
  cNames[idx] <- paste0("DS:", cNames[idx])
  
  colnames(dataMatrix) <- cNames
  
  for(ii in 1:nrow(dataMatrix)){
    
    cc <- conditions[[ii]]
    
    for(kk in 1:length(sites)){
      
      nn <- dataGMM@IDmap$dataID[which(dataGMM@IDmap$S.cc==sites[kk])]
      idxMeas <- which(names(dataGMM@res)%in%nn)
      
      if(length(idxMeas) > 0){
        
        score <- min(as.numeric(dataGMM@res[idxMeas][[1]][which(rownames(dataGMM@res[idxMeas][[1]])%in%cc), 1]))
        
        dataMatrix[ii, which(cNames==paste0("DS:", sites[kk]))] <- score
        
      }
      
    }
    
  }
  
  if(length(which(is.na(dataMatrix)))>0){dataMatrix[which(is.na(dataMatrix))] <- 0}
  if(length(intersect(which(dataMatrix < 0), which(is.infinite(dataMatrix))))>0){dataMatrix[intersect(which(dataMatrix < 0), which(is.infinite(dataMatrix)))] <- -9999}
  if(length(dataMatrix[intersect(which(dataMatrix > 0), which(is.infinite(dataMatrix)))] <- 9999)>0){dataMatrix[intersect(which(dataMatrix > 0), which(is.infinite(dataMatrix)))] <- 9999}
  
  for(ii in 1:length(targets)){
    
    idxTG[[length(idxTG)+1]] <- which(cNames%in%paste0("DN:", targets[[ii]]))
    idxDN[[length(idxDN)+1]] <- setdiff(1:length(cNames), unique(c(idxDS[[ii]], idxTG[[ii]])))
    
  }
  
  res <- list(dataMatrix=dataMatrix, tgID=idxTG, dnID=idxDN, dsID=idxDS, species=unique(c(TG, setdiff(pknList@species, c(TG, sites)), sites)))
  
  return(res)
  
}