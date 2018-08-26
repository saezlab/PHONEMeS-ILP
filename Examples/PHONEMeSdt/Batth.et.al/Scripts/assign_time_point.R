#load the packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(CellNOptR)

load(file = "resList.RData")

load(file = "dataGMM.RData")

measuredSites <- GMM.ID$S.cc

interactions <- matrix(data = , nrow = 1, ncol = 4)
colnames(interactions) <- c("Source", "Weight", "Target", "tp")

##
for(ii in 1:length(resList)){
  
  for(jj in 1:nrow(resList[[ii]][[1]])){
    
    toBind <- t(as.matrix(c(resList[[ii]][[1]][jj, 1], "1", resList[[ii]][[1]][jj, 3], "3min")))
    colnames(toBind) <- colnames(interactions)
    
    interactions <- rbind(interactions, toBind)
    
  }
  
}

interactions <- interactions[-1, ]
interactions <- unique(interactions)

for(ii in 1:nrow(interactions)){
  
  cnt <- 0
  for(jj in 1:length(resList)){
    
    idx1 <- which(resList[[jj]][[1]][, 1]==interactions[ii, 1])
    idx2 <- which(resList[[jj]][[1]][, 3]==interactions[ii, 3])
    
    idx <- intersect(x = idx1, y = idx2)
    
    if(length(idx) > 0){
      
      cnt <- cnt + 1
      
    }
    
  }
  
  interactions[ii, 2] <- as.character(cnt)
  
}

##
for(ii in 1:length(resList)){
  
  for(jj in 1:nrow(resList[[ii]][[2]])){
    
    idx1 <- which(interactions[, 1]==resList[[ii]][[2]][jj, 1])
    idx2 <- which(interactions[, 3]==resList[[ii]][[2]][jj, 3])
    
    idx <- intersect(x = idx1, y = idx2)
    
    if(length(idx)==0){
      
      toBind <- t(as.matrix(c(resList[[ii]][[2]][jj, 1], "1", resList[[ii]][[2]][jj, 3], "15min")))
      colnames(toBind) <- colnames(interactions)
      
      interactions <- rbind(interactions, toBind)
      
    }
    
  }
  
}

interactions <- unique(interactions)

tpIDX <- which(interactions[, 4]=="15min")

for(ii in 1:length(tpIDX)){
  
  cnt <- 0
  for(jj in 1:length(resList)){
    
    idx1 <- which(resList[[jj]][[2]][, 1]==interactions[tpIDX[ii], 1])
    idx2 <- which(resList[[jj]][[2]][, 3]==interactions[tpIDX[ii], 3])
    
    idx <- intersect(x = idx1, y = idx2)
    
    if(length(idx) > 0){
      
      cnt <- cnt + 1
      
    }
    
  }
  
  interactions[tpIDX[ii], 2] <- as.character(cnt)
  
}

##
# interactions <- interactions[-which(interactions[, 1]%in%c("CAMK2D", "CAMK2D_T331")), ]
write.table(x = interactions, file = "interactions-all.txt", quote = FALSE, sep = "\t", row.names = FALSE)

threshWeight <- 40
measuredSpecies <- GMM.ID$S.cc

gg <- graph_from_data_frame(d = as.data.frame(interactions[, c(1, 3)]), directed = TRUE)
adj <- get.adjacency(graph = gg)

paths2add <- list()

for(ii in 1:length(measuredSpecies)){
  
  if(measuredSpecies[ii] %in% rownames(adj)){
    
    path <- all_simple_paths(graph = gg, from = which(rownames(adj)=="PDGFR"), to = which(rownames(adj)==measuredSpecies[ii]))
    
    for(i in 1:length(path)){
      
      cnt <- 0
      for(j in 1:(length(path[[i]])-1)){
        
        ss <- rownames(adj)[path[[i]][j]]
        tt <- rownames(adj)[path[[i]][j+1]]
        
        idx <- intersect(x = which(interactions[, 1]==ss), y = which(interactions[, 3]==tt))
        
        if(as.numeric(interactions[idx, 2]) < threshWeight){
          
          cnt <- cnt + 1
          
        }
        
      }
      
      if(cnt==0){
        
        paths2add[[length(paths2add)+1]] <- path[[i]]
        
      }
      
    }
    
  }
  
}

interactionMatrix <- matrix(data = , nrow = 1, ncol = 4)
colnames(interactionMatrix) <- c("Source", "Weight", "Target", "tp")

for(ii in 1:length(paths2add)){
  
  for(jj in 1:(length(paths2add[[ii]])-1)){
    
    ss <- rownames(adj)[paths2add[[ii]][jj]]
    tt <- rownames(adj)[paths2add[[ii]][jj+1]]
    
    idx <- intersect(x = which(interactions[, 1]==ss), y = which(interactions[, 3]==tt))
    
    cc <- t(interactions[idx, ])
    
    interactionMatrix <- unique(rbind(interactionMatrix, cc))
    
  }
  
}

interactionMatrix <- interactionMatrix[-1, ]

write.table(x = interactionMatrix, file = "interactions-cutoff.txt", quote = FALSE, sep = "\t", row.names = FALSE)

##
nodesAttributes <- matrix(data = , nrow = length(GMM.ID$S.cc)+1, ncol = 2)
colnames(nodesAttributes) <- c("Sites", "nodesP")

nodesAttributes[1, 1] <- "PDGFR"
nodesAttributes[1, 2] <- "D"

nodesAttributes[2:nrow(nodesAttributes), 1] <- GMM.ID$S.cc
nodesAttributes[2:nrow(nodesAttributes), 2] <- "P"

write.table(x = nodesAttributes, file = "nodesAttributes.txt", quote = FALSE, sep = "\t", row.names = FALSE)
