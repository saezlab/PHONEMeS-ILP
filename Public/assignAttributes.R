assignAttributes <- function(sif = sif, dataGMM = dataGMM, targets = targets, writeAttr = TRUE, fileName = "nodesAttributes.txt"){
  
  species <- unique(c(sif[, 1], sif[, 3]))
  
  nodesAttributes <- matrix(data = , nrow = length(species), ncol = 2)
  colnames(nodesAttributes) <- c("Species", "nodesP")
  
  nodesAttributes[, 1] <- species
  nodesAttributes[, 2] <- ""
  
  tt <- unique(unlist(targets))
  nodesAttributes[which(nodesAttributes[, 1]%in%tt), 2] <- "D"
  
  ss <- dataGMM@IDmap$S.cc
  nodesAttributes[which(nodesAttributes[, 1]%in%ss), 2] <- "P"
  
  if(writeAttr){
    
    write.table(x = nodesAttributes, file = fileName, quote = FALSE, sep = "\t", row.names = FALSE)
    
  }
  
  return(nodesAttributes)
  
}