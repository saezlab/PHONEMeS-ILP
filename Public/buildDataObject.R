buildDataObject <- function(ttList = ttList, pThresh = 0.1, pValIdx = 6, measIdx = 1, fcIdx = 2, organism = "HUMAN"){
  
  print("Building inputs for PHONEMeS...")
  
  tt.p <- matrix(data = , nrow = length(ttList), ncol = nrow(ttList[[1]]))
  tt.lo <- matrix(data = , nrow = length(ttList), ncol = nrow(ttList[[1]]))
  tt.c <- matrix(data = , nrow = length(ttList), ncol = nrow(ttList[[1]]))
  tt.s <- matrix(data = , nrow = length(ttList), ncol = nrow(ttList[[1]]))
  tt.fc <- matrix(data = , nrow = length(ttList), ncol = nrow(ttList[[1]]))
  
  rownames(tt.p) <- names(ttList)
  rownames(tt.lo) <- names(ttList)
  rownames(tt.c) <- names(ttList)
  rownames(tt.s) <- names(ttList)
  rownames(tt.fc) <- names(ttList)
  
  colnames(tt.p) <- ttList[[1]][, measIdx]
  colnames(tt.lo) <- ttList[[1]][, measIdx]
  colnames(tt.c) <- ttList[[1]][, measIdx]
  colnames(tt.s) <- ttList[[1]][, measIdx]
  colnames(tt.fc) <- ttList[[1]][, measIdx]
  
  for(ii in 1:length(ttList)){
    
    for(jj in 1:nrow(ttList[[ii]])){
      
      if((is.na(ttList[[ii]][jj, pValIdx])) || (is.na(ttList[[ii]][jj, measIdx])) || (is.na(ttList[[ii]][jj, fcIdx]))){
        
        tt.p[ii, which(colnames(tt.p)==ttList[[ii]][jj, measIdx])] <- "1"
        tt.lo[ii, which(colnames(tt.p)==ttList[[ii]][jj, measIdx])] <- "100"
        tt.c[ii, which(colnames(tt.p)==ttList[[ii]][jj, measIdx])] <- "C"
        tt.s[ii, jj] <- "OK"
        tt.fc[ii, which(colnames(tt.p)==ttList[[ii]][jj, measIdx])] <- "0"
        
      } else {
        
        tt.p[ii, which(colnames(tt.p)==ttList[[ii]][jj, measIdx])] <- ttList[[ii]][jj, pValIdx]
        tt.lo[ii, which(colnames(tt.p)==ttList[[ii]][jj, measIdx])] <- log(x = (ttList[[ii]][jj, pValIdx])/pThresh)
        if(ttList[[ii]][jj, pValIdx] < pThresh){
          tt.c[ii, which(colnames(tt.p)==ttList[[ii]][jj, measIdx])] <- "P"
        } else {
          tt.c[ii, which(colnames(tt.p)==ttList[[ii]][jj, measIdx])] <- "C"
        }
        tt.s[ii, jj] <- "OK"
        tt.fc[ii, which(colnames(tt.p)==ttList[[ii]][jj, measIdx])] <- ttList[[ii]][jj, fcIdx]
        
      }
      
    }
    
  }
  
  ##
  GMM<-vector("list", length = nrow(ttList[[1]]))
  names(GMM)<-ttList[[1]][, measIdx]
  
  GMM.wFC<-vector("list", length = nrow(ttList[[1]]))
  names(GMM.wFC)<-ttList[[1]][, measIdx]
  
  for(i in 1:length(GMM)){
    
    GMM[[i]]<-matrix(data = , nrow = length(ttList), ncol = 4)
    GMM.wFC[[i]]<-matrix(data = , nrow = length(ttList), ncol = 5)
    
    for(j in 1:length(ttList)){
      
      GMM[[i]][j, 1] <- as.character(tt.lo[j, names(GMM)[i]])
      GMM[[i]][j, 2] <- as.character(tt.c[j, names(GMM)[i]])
      GMM[[i]][j, 3] <- as.character(tt.p[j, names(GMM)[i]])
      GMM[[i]][j, 4] <- as.character(tt.s[j, names(GMM)[i]])
      
      GMM.wFC[[i]][j, 1] <- as.character(tt.lo[j, names(GMM.wFC)[i]])
      GMM.wFC[[i]][j, 2] <- as.character(tt.c[j, names(GMM.wFC)[i]])
      GMM.wFC[[i]][j, 3] <- as.character(tt.p[j, names(GMM.wFC)[i]])
      GMM.wFC[[i]][j, 4] <- as.character(tt.s[j, names(GMM.wFC)[i]])
      GMM.wFC[[i]][j, 5] <- as.character(tt.fc[j, names(GMM.wFC)[i]])
      
    }
    
    colnames(GMM[[i]])<-c("Indiv", "clus","FCvCaPval","status")
    rownames(GMM[[i]])<-names(ttList)
    
    colnames(GMM.wFC[[i]])<-c("Indiv", "clus","FCvCaPval","status", "logFC")
    rownames(GMM.wFC[[i]])<-names(ttList)
    
  }
  
  ##
  data.IDmap <- matrix(data = , nrow = nrow(ttList[[1]]), ncol = 3)
  colnames(data.IDmap) <- c("dataID", "UPID", "S.cc")
  data.IDmap[, 1] <- ttList[[1]][, measIdx]
  idx <- which(grepl(pattern = "_HUMAN", x = ttList[[1]][, measIdx]))
  if(length(idx)==0){
    
    data.IDmap[, 2] <- paste0(sapply(strsplit(x = ttList[[1]][, measIdx], split = "_", fixed = TRUE), "[[", 1), "_", organism)
    data.IDmap[, 3] <- paste0(sapply(strsplit(x = ttList[[1]][, measIdx], split = "_", fixed = TRUE), "[[", 1), "_", organism, "_", sapply(strsplit(x = ttList[[1]][, measIdx], split = "_", fixed = TRUE), "[[", 2))
    
  }
  
  GMM.ID <- as.data.frame(data.IDmap)
  
  save(list=c("GMM.ID", "GMM","GMM.wFC"), file="dataGMM.RData")
  
  print("DONE!!!")
  
  return(list=c("GMM.ID", "GMM","GMM.wFC"))
  
}