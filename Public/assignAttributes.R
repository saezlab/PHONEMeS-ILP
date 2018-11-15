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

# Assigning node attributes:
#   P <- measured site (red hexagonal shape)
#   D <- target node (green rombe shape)

assignAttributes <- function(sif = sif, dataGMM = dataGMM, targets = targets, writeAttr = TRUE){
  
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
    
    write.table(x = nodesAttributes, file = "nodesAttributes.txt", quote = FALSE, sep = "\t", row.names = FALSE)
    
  }
  
  return(nodesAttributes)
  
}