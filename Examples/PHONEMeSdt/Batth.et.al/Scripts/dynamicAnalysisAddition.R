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
# 13:00 23/03/2018
# This script is used to assign the attributes for the newly added edges after
# each time-point.

library(readr)
timePoints <- 6

for(ii in 1:timePoints){
  
  if(ii==1){
    
    currSIF <- read_delim(paste0("resultSIF_tp", ii, ".txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
    allSIF <- currSIF
    
  }
  else{
    
    currSIF <- read_delim(paste0("resultSIF_tp", ii, ".txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
    idx <- c()
    for(i in 1:nrow(currSIF)){
      
      if(nrow(merge(currSIF[i, ],allSIF))>0){
        
        idx <- c(idx, i)
        
      }
      
    }
    
    if(length(idx) < nrow(currSIF)){
      
      diffSIF <- currSIF[-idx, ]
      write.table(diffSIF, file = paste0("diffSIF_tp_", ii-1, "-", ii, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
      
    }
    else{
      
      write.table(matrix(, nrow = 1, ncol = 3), file = paste0("diffSIF_tp_", ii-1, "-", ii, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
      
    }
    
    allSIF <- unique(rbind(allSIF, currSIF))
    
  }
  
}

#######
# Assigning Edge attributes

#
resultSIF_tp1 <- read_delim("resultSIF_tp1.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
resultSIF_tp1[, 2] <- "100"
colnames(resultSIF_tp1) <- c("Source", "Width", "Target")
write.table(resultSIF_tp1, file = "resultSIF_tp1_added.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#
resultSIF_tp2 <- read_delim("resultSIF_tp2.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
resultSIF_tp2[, 2] <- "10"
resultSIF_tp2 <- as.matrix(resultSIF_tp2)
diffSIF_tp_1_2 <- read_delim("diffSIF_tp_1-2.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
diffSIF_tp_1_2 <- as.matrix(diffSIF_tp_1_2)
idx <- c()
for(i in 1:nrow(diffSIF_tp_1_2)){
  
  kk <- intersect(which(resultSIF_tp2[, 1]==diffSIF_tp_1_2[i, 1]), which(resultSIF_tp2[, 3]==diffSIF_tp_1_2[i, 3]))
  
  idx <- c(idx, kk)
  
}
resultSIF_tp2[idx, 2] <- "100"
colnames(resultSIF_tp2) <- c("Source", "Width", "Target")
write.table(resultSIF_tp2, file = "resultSIF_tp2_added.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# #
# resultSIF_tp3 <- read_delim("resultSIF_tp3.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# resultSIF_tp3[, 2] <- "10"
# resultSIF_tp3 <- as.matrix(resultSIF_tp3)
# diffSIF_tp_2_3 <- read_delim("diffSIF_tp_2-3.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# diffSIF_tp_2_3 <- as.matrix(diffSIF_tp_2_3)
# idx <- c()
# for(i in 1:nrow(diffSIF_tp_2_3)){
#   
#   kk <- intersect(which(resultSIF_tp3[, 1]==diffSIF_tp_2_3[i, 1]), which(resultSIF_tp3[, 3]==diffSIF_tp_2_3[i, 3]))
#   
#   idx <- c(idx, kk)
#   
# }
# resultSIF_tp3[idx, 2] <- "100"
# colnames(resultSIF_tp3) <- c("Source", "Width", "Target")
# write.table(resultSIF_tp3, file = "resultSIF_tp3_added.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# 
# #
# resultSIF_tp4 <- read_delim("resultSIF_tp4.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# resultSIF_tp4[, 2] <- "10"
# resultSIF_tp4 <- as.matrix(resultSIF_tp4)
# diffSIF_tp_3_4 <- read_delim("diffSIF_tp_3-4.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# diffSIF_tp_3_4 <- as.matrix(diffSIF_tp_3_4)
# idx <- c()
# for(i in 1:nrow(diffSIF_tp_3_4)){
#   
#   kk <- intersect(which(resultSIF_tp4[, 1]==diffSIF_tp_3_4[i, 1]), which(resultSIF_tp4[, 3]==diffSIF_tp_3_4[i, 3]))
#   
#   idx <- c(idx, kk)
#   
# }
# resultSIF_tp4[idx, 2] <- "100"
# colnames(resultSIF_tp4) <- c("Source", "Width", "Target")
# write.table(resultSIF_tp4, file = "resultSIF_tp4_added.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# 
# #
# resultSIF_tp5 <- read_delim("resultSIF_tp5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# resultSIF_tp5[, 2] <- "10"
# resultSIF_tp5 <- as.matrix(resultSIF_tp5)
# diffSIF_tp_4_5 <- read_delim("diffSIF_tp_4-5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# diffSIF_tp_4_5 <- as.matrix(diffSIF_tp_4_5)
# idx <- c()
# for(i in 1:nrow(diffSIF_tp_4_5)){
#   
#   kk <- intersect(which(resultSIF_tp5[, 1]==diffSIF_tp_4_5[i, 1]), which(resultSIF_tp5[, 3]==diffSIF_tp_4_5[i, 3]))
#   
#   idx <- c(idx, kk)
#   
# }
# resultSIF_tp5[idx, 2] <- "100"
# colnames(resultSIF_tp5) <- c("Source", "Width", "Target")
# write.table(resultSIF_tp5, file = "resultSIF_tp5_added.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# 
# #
# resultSIF_tp6 <- read_delim("resultSIF_tp6.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# resultSIF_tp6[, 2] <- "10"
# resultSIF_tp6 <- as.matrix(resultSIF_tp6)
# diffSIF_tp_5_6 <- read_delim("diffSIF_tp_5-6.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# diffSIF_tp_5_6 <- as.matrix(diffSIF_tp_5_6)
# idx <- c()
# for(i in 1:nrow(diffSIF_tp_5_6)){
#   
#   kk <- intersect(which(resultSIF_tp6[, 1]==diffSIF_tp_5_6[i, 1]), which(resultSIF_tp6[, 3]==diffSIF_tp_5_6[i, 3]))
#   
#   idx <- c(idx, kk)
#   
# }
# resultSIF_tp6[idx, 2] <- "100"
# colnames(resultSIF_tp6) <- c("Source", "Width", "Target")
# write.table(resultSIF_tp6, file = "resultSIF_tp6_added.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# 
# 
# 
# 
