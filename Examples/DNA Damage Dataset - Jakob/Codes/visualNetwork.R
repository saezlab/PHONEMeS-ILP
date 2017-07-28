measuredSites <- c()
for(i in 1:length(dataMatrix$dsID)){
  
  if(dataMatrix$dataMatrix[1, dataMatrix$dsID[i]] < 0){
    
    measuredSites <- c(measuredSites, strsplit(colnames(dataMatrix$dataMatrix)[dataMatrix$dsID[i]], split = ":")[[1]][2])
    
  }
  
}


allNodes <- unique(c(as.character(resultsSIF[, 1]), as.character(resultsSIF[, 3])))

allNodesMatrix <- matrix(, nrow = length(measuredSites)+1, ncol = 2)
allNodesMatrix[1, 1] <- "ATM_HUMAN"
allNodesMatrix[1, 2] <- "D"
for(i in 1:length(measuredSites)){
  
  allNodesMatrix[i+1, 1] <- measuredSites[i]
  
  allNodesMatrix[i+1, 2] <- "P"
  
}

write.table(x = allNodesMatrix, file = "nodesAttributes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
