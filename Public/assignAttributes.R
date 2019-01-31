assignAttributes <- function(sif = sif, dataGMM = dataGMM, targets = targets, writeAttr = TRUE){
  
  species <- unique(c(sif[, 1], sif[, 3]))
  
  nodesAttributes <- matrix(data = , nrow = length(species), ncol = 3)
  colnames(nodesAttributes) <- c("Species", "nodesP", "nodesInf")
  
  nodesAttributes[, 1] <- species
  nodesAttributes[, 2] <- ""
  
  tt <- unique(unlist(targets))
  nodesAttributes[which(nodesAttributes[, 1]%in%tt), 2] <- "D"
  
  ss <- dataGMM@IDmap$S.cc
  nodesAttributes[which(nodesAttributes[, 1]%in%ss), 2] <- "P"
  
  sites <- nodesAttributes[which(grepl(pattern = ".", x = nodesAttributes[, 1], fixed = TRUE)), 1]
  sitesInferred <- setdiff(x = sites, y = nodesAttributes[which(nodesAttributes[, 2]=="P"), 1])
  
  nodesAttributes[, 3] <- "Normal"
  nodesAttributes[which(nodesAttributes[, 1]%in%sitesInferred), 3] <- "Small"
  
  protein <- c()
  offset <- c()
  
  for(ii in 1:length(sitesInferred)){
    
    protein <- c(protein, strsplit(x = sitesInferred[ii], split = "_HUMAN.")[[1]][1])
    offset <- c(offset, strsplit(x = sitesInferred[ii], split = "_HUMAN.")[[1]][2])
    
  }
  
  mapping <- matrix(data = , nrow = length(sitesInferred), ncol = 2)
  mapping[, 1] <- sitesInferred
  
  for(ii in 1:nrow(mapping)){
    
    pp <- strsplit(x = mapping[ii, 1], split = "_", fixed = TRUE)[[1]][1]
    idx <- which(protein==pp)
    
    for(jj in 1:length(idx)){
      
      if(jj==1){
        
        mapping[ii, 2] <- offset[idx[jj]]
        
      } else {
        
        mapping[ii, 2] <- paste0(mapping[ii, 2], ";", offset[idx[jj]])
        
      }
      
    }
    
  }
  
  returnAttributes <- nodesAttributes
  for(ii in 1:nrow(returnAttributes)){
    
    if(returnAttributes[ii, 1]%in%mapping[, 1]){ returnAttributes[ii, 1] <- mapping[which(mapping[, 1]==returnAttributes[ii, 1]), 2]}
    
  }
  
  returnAttributes <- unique(returnAttributes)
  
  returnSIF <- sif
  for(ii in 1:nrow(returnSIF)){
    
    if(returnSIF[ii, 1]%in%mapping[, 1]){returnSIF[ii, 1] <- mapping[which(mapping[, 1]==returnSIF[ii, 1]), 2]}
    if(returnSIF[ii, 3]%in%mapping[, 1]){returnSIF[ii, 3] <- mapping[which(mapping[, 1]==returnSIF[ii, 3]), 2]}
    
  }
  
  temp <- returnSIF
  
  returnSIF[, 2] <- 0
  returnSIF <- unique(returnSIF)
  for(ii in 1:nrow(returnSIF)){
    
    idx1 <- which(temp[, 1]==returnSIF[ii, 1])
    idx2 <- which(temp[, 3]==returnSIF[ii, 3])
    
    idx <- intersect(x = idx1, y = idx2)
    
    returnSIF[ii, 2] <- sum(temp[idx, 2])
    
  }
  
  colnames(returnSIF) <- c("Source", "f50", "Target")
  
  if(writeAttr){
    
    write.table(x = nodesAttributes, file = "nodesAttributes.txt", quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(x = returnSIF, file = "returnSIF.txt", quote = FALSE, sep = "\t", row.names = FALSE)
    
  }
  
  resList <- list()
  resList[[1]] <- returnSIF
  resList[[2]] <- nodesAttributes
  
  return(resList)
  
}