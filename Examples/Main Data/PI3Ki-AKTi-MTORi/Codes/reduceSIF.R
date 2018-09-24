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

# Function to reduce the the final PHONEMeS network and select only some specific
# most likely mechanisms for visualization

reduceSIF <- function(sif = sif, targets = targets, dataGMM = dataGMM, cutoff = 0.1){
  
  targets <- unique(unlist(targets))
  allSpecies <- intersect(x = allSpecies <- unique(c(sif[, 1], sif[, 3])), dataGMM@IDmap$S.cc)
  
  maxCnt <- max(as.numeric(sif[, 2]))
  cntCutoff <- maxCnt*cutoff
  
  gg <- graph_from_data_frame(d = as.data.frame(x = sif[, c(1, 3)]), directed = TRUE)
  adj <- get.adjacency(graph = gg)
  
  returnSIF <- matrix(data = , nrow = 1, ncol = 3)
  colnames(returnSIF) <- colnames(sif)
  for(ii in 1:length(targets)){
    
    for(jj in 1:length(allSpecies)){
      
      sP <- get.all.shortest.paths(graph = gg, from = which(rownames(adj)==targets[ii]), to = which(rownames(adj)==allSpecies[jj]))
      
      if(length(sP[[1]]) > 0){
        
        for(kk in 1:length(sP[[1]])){
          
          cnt <- 0
          temp <- matrix(data = , nrow = 1, ncol = 3)
          colnames(temp) <- colnames(returnSIF)
          for(ll in 1:(length(sP[[1]][[kk]])-1)){
            
            idx1 <- which(sif[, 1]==rownames(adj)[sP[[1]][[kk]][ll]])
            idx2 <- which(sif[, 3]==rownames(adj)[sP[[1]][[kk]][ll+1]])
            
            idx <- intersect(x = idx1, y = idx2)
            
            if(as.numeric(sif[idx, 2]) < cntCutoff){
              
              cnt <- cnt + 1
              
            }
            
            temp <- rbind(temp, sif[idx, ])
            
          }
          
          if(cnt == 0){
            
            returnSIF <- rbind(returnSIF, temp[-1, ])
            
          }
          
        }
        
      }
      
    }
    
  }
  
  returnSIF <- returnSIF[-1, ]
  returnSIF <- unique(returnSIF)
  
  return(returnSIF)
  
}