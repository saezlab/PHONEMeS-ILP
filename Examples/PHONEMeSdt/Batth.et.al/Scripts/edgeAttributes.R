load("dataGMM.RData")

measuredSites <- c()
for(i in 1:length(names(GMM))){
  
  measuredSites <- c(measuredSites, as.character(GMM.ID$S.cc)[which(as.character(GMM.ID$dataID)==names(GMM)[i])])
  
}

for(ii in 1:5){
  
  currSIF <-  read_delim(paste0("resultsSIF_tp", ii, ".txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
  mat <- matrix(, nrow = nrow(currSIF), ncol = 3)
  mat[, 1] <- currSIF$V1
  mat[, 3] <- currSIF$V3
  
  if(ii==1){
    
    mat[, 2] <- rep("yes", nrow(mat))
    write.table(x = mat, file = "resultsSIF_tp1_edges.txt", quote = FALSE, row.names = FALSE, sep = "\t")
    
  }
  else{
    
    mat[, 2] <- rep("yes", nrow(mat))
    for(i in 1:nrow(mat)){
      
      if(mat[i, 1]%in%measuredSites){
        
        idx <- which(names(GMM)==as.character(GMM.ID$dataID)[which(as.character(GMM.ID$S.cc)==mat[i, 1])])
        
        if(GMM[[idx]][ii, 2]=="C"){
          
          mat[i, 2] <- "no"
          
        }
        
      }
      
      if(mat[i, 3]%in%measuredSites){
        
        idx <- which(names(GMM)==as.character(GMM.ID$dataID)[which(as.character(GMM.ID$S.cc)==mat[i, 3])])
        
        if(GMM[[idx]][ii, 2]=="C"){
          
          mat[i, 2] <- "no"
          
        }
        
      }
      
    }
    
    write.table(x = mat, file = paste0("resultsSIF_tp", ii, "_edges.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
    
  }
  
}