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
# 13:24 23/03/2018
# This script is used to assign and combine all the edge attributes

for(ii in 1:2){
  
  aSIF <- read_delim(paste0("resultSIF_tp", ii, "_added.txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
  rSIF <- read_delim(paste0("resultsSIF_removed_tp", ii, ".txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
  
  sif <- matrix(data = , nrow = nrow(aSIF), ncol = 4)
  
  sif[, 1] <- as.character(as.matrix(aSIF[, 1]))
  sif[, 2] <- as.character(as.matrix(rSIF[, 3]))
  
  sif[, 3] <- as.character(as.matrix(aSIF[, 2]))
  sif[, 4] <- as.character(as.matrix(rSIF[, 2]))
  
  colnames(sif) <- c("Source", "Target", "Width", "State")
  
  idx <- intersect(which(sif[, 3]=="100"), which(sif[, 4]=="W"))
  if(length(idx) > 0){
    
    for(i in 1:length(idx)){sif[idx[i], 3] <- "10"}
    
  }
  
  write.table(x = sif, file = paste0("resultSIF_wAttrib_tp", ii, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
}