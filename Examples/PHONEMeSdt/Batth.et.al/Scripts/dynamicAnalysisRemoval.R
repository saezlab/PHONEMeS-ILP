#  Copyright (c) 2018 - RWTH Aachen University
#
#  File author(s): Enio Gjerga
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  email: enio.gjerga@gmai.com
#
##############################################################################
# 13:23 23/03/2018
# This script is used to assign the attributes for the removed pathways after
# each time-point.

currSIF <- read_delim(paste0("resultSIF_tp", 1, ".txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
currSIF[, 2] <- "S"
toSave <- currSIF

colnames(toSave) <- c("Source", "Interaction", "Target")

write.table(x = toSave, file = paste0("resultsSIF_removed_tp1.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

for(ii in 2:2){
  
  currSIF <- read_delim(paste0("resultSIF_tp", ii, ".txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
  
  species <- unique(as.character(unlist(c(currSIF[, 1], currSIF[, 3]))))
  
  sites <- c()
  for(i in 1:length(species)){
    
    if(length(strsplit(x = species[i], split = "_", fixed = TRUE)[[1]])==2){
      
      sites <- c(sites, species[i])
      
    }
    
  }
  
  sites <- intersect(sites, GMM.ID$S.cc)
  
  nameSites <- c()
  for(i in 1:length(sites)){
    
    nameSites <- c(nameSites, as.character(GMM.ID$dataID[which(GMM.ID$S.cc==sites[i])]))
    
  }
  
  perturbations <- c()
  for(i in 1:length(nameSites)){
    
    idx <- which(names(GMM.wFC)==nameSites[i])
    if(GMM.wFC[[idx]][ii, 2]=="P"){
      
      print(GMM.wFC[[idx]][ii, 2])
      perturbations <- c(perturbations, as.character(GMM.ID$S.cc[which(GMM.ID$dataID==nameSites[i])]))
      
    }
    
  }
  
  controls <- c()
  for(i in 1:length(nameSites)){
    
    idx <- which(names(GMM.wFC)==nameSites[i])
    if(GMM.wFC[[idx]][ii, 2]=="C"){
      
      controls <- c(controls, as.character(GMM.ID$S.cc[which(GMM.ID$dataID==nameSites[i])]))
      
    }
    
  }
  
  controls <- setdiff(controls, perturbations)
  
  prevSIF <- read_delim(paste0("resultSIF_tp", ii-1, ".txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # Find ctrl paths
  
  g1 <- graph_from_data_frame(d = currSIF[, c(1, 3)], directed = TRUE)
  
  adj1 <- get.adjacency(graph = g1)
  
  pathsCtrl <- list()
  
  for(i in 1:length(controls)){
    
    # pathsCtrl[[length(pathsCtrl)+1]] <- get.all.shortest.paths(graph = g1, from = which(rownames(adj1)=="EDNRB_HUMAN"), to = which(rownames(adj1)==controls[i]))
    pathsCtrl[[length(pathsCtrl)+1]] <- all_simple_paths(graph = g1, from = which(rownames(adj1)=="EDNRB_HUMAN"), to = which(rownames(adj1)==controls[i]))
    
  }
  
  pathsPert <- list()
  
  for(i in 1:length(perturbations)){
    
    # pathsPert[[length(pathsPert)+1]] <- get.all.shortest.paths(graph = g1, from = which(rownames(adj1)=="EDNRB_HUMAN"), to = which(rownames(adj1)==perturbations[i]))
    pathsPert[[length(pathsPert)+1]] <- all_simple_paths(graph = g1, from = which(rownames(adj1)=="EDNRB_HUMAN"), to = which(rownames(adj1)==perturbations[i]))
    
  }
  
  ctrlSIF <- matrix(data = , nrow = 1, ncol = 3)
  
  for(i in 1:length(pathsCtrl)){
    
    if(length(pathsCtrl[[i]]) > 0){
      
      for(j in 1:length(pathsCtrl[[i]])){
        
        if(length(pathsCtrl[[i]][[j]]) > 1){
          
          pp <- pathsCtrl[[i]][[j]]
          
          for(k in 1:(length(pp)-1)){
            
            ctrlSIF <- rbind(ctrlSIF, c(rownames(adj1)[pp[k]], "1", rownames(adj1)[pp[k+1]]))
            
          }
          
        }
        
      }
      
    }
    
    # if(length(pathsCtrl[[i]][[i]]) > 0){
    #   
    #   pp <- pathsCtrl[[i]][[1]]
    #   
    #   for(j in 1:(length(pp)-1)){
    #     
    #     ctrlSIF <- rbind(ctrlSIF, c(rownames(adj1)[pp[j]], "1", rownames(adj1)[pp[j+1]]))
    #     
    #   }
    #   
    # }
    
  }
  ctrlSIF <- unique(ctrlSIF[-1, ])
  
  pertSIF <- matrix(data = , nrow = 1, ncol = 3)
  for(i in 1:length(pathsPert)){
    
    if(length(pathsPert[[i]]) > 0){
      
      for(j in 1:length(pathsPert[[i]])){
        
        if(length(pathsPert[[i]][[j]]) > 1){
          
          pp <- pathsPert[[i]][[j]]
          
          for(k in 1:(length(pp)-1)){
            
            pertSIF <- rbind(pertSIF, c(rownames(adj1)[pp[k]], "1", rownames(adj1)[pp[k+1]]))
            
          }
          
        }
        
      }
      
    }
    
    # if(length(pathsPert[[i]][[1]]) > 0){
    #   
    #   pp <- pathsPert[[i]][[1]][[1]]
    #   
    #   for(j in 1:(length(pp)-1)){
    #     
    #     pertSIF <- rbind(pertSIF, c(rownames(adj1)[pp[j]], "1", rownames(adj1)[pp[j+1]]))
    #     
    #   }
    #   
    # }
    
  }
  pertSIF <- unique(pertSIF[-1, ])
  
  diffSIF <- matrix(data = , nrow = 1, ncol = 3)
  
  idx <- c()
  for(i in 1:nrow(ctrlSIF)){
    
    cnt <- 0
    for(j in 1:nrow(pertSIF)){
      
      if(ctrlSIF[i, 1]==pertSIF[j, 1] && ctrlSIF[i, 3]==pertSIF[j, 3]){
        
        cnt <- cnt + 1
        
      }
      
    }
    
    if(cnt==0){
      
      idx <- c(idx, i)
      
    }
    
  }
  
  diffSIF <- rbind(diffSIF, ctrlSIF[idx, ])[-1, ]
  
  toSave <- currSIF
  toSave[, 2] <- rep("S")
  
  for(i in 1:nrow(toSave)){
    
    cnt <- 0
    for(j in 1:nrow(diffSIF)){
      
      if(toSave[i, 1]==diffSIF[j, 1] && toSave[i, 3]==diffSIF[j, 3]){
        
        toSave[i, 2] <- "W"
        
      }
      
    }
    
  }
  
  colnames(toSave) <- c("Source", "Interaction", "Target")
  
  write.table(x = toSave, file = paste0("resultsSIF_removed_tp", ii, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
}
