buildInputs <- function(tableTopList = tableTopList, organism = organism, pThresh = pThres, mappingTable = mappingTable){
  
  ###
  # Build GMM.ID object
  
  if(is.null(mappingTable)){
    
    print("Preparing Mapping Table..")
    if(is.null(organism)){
      
      uniprot <- read_delim("uniprot-all-HUMAN.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
      
    } else {
      
      uniprot <- read_delim(paste0("uniprot-all-", organism, ".tab"), "\t", escape_double = FALSE, trim_ws = TRUE)
      
    }
    
    species <- c()
    for(i in 1:length(tableTopList)){
      
      species <- c(species, tableTopList[[i]]$X1)
      
    }
    
    species <- unique(species)
    
    proteins <- unlist(lapply(strsplit(x = species, split = "_", fixed = TRUE), '[[', 1))
    residues <- unlist(lapply(strsplit(x = species, split = "_", fixed = TRUE), '[[', 2))
    
    proteins <- proteins[-which(residues=="")]
    residues <- residues[-which(residues=="")]
    
    Table.ID <- matrix(data = , nrow = length(proteins), ncol = 6)
    colnames(Table.ID) <- c("dataID", "UPID", "site", "pos", "res", "S.cc")
    
    Table.ID[, 1] <- paste0(proteins, "_", residues)
    Table.ID[, 3] <- residues
    for(i in 1:nrow(Table.ID)){
      
      Table.ID[i, 2] <- as.character(uniprot$`Entry name`[which(uniprot$Entry==proteins[i])[1]])
      Table.ID[i, 4] <- substr(x = residues[i], start = 2, stop = nchar(residues[i]))
      Table.ID[i, 5] <- substr(x = residues[i], start = 1, stop = 1)
      Table.ID[i, 6] <- paste0(Table.ID[i, 2], "_", residues[i])
      
    }
    
    Table.ID <- Table.ID[complete.cases(Table.ID), ]
    
    GMM.ID <- as.data.frame(Table.ID)
    
  } else {
    
    GMM.ID <- as.data.frame(mappingTable)
    Table.ID <- as.matrix(mappingTable)
    
  }
  
  ###
  # Build data object
  
  print("Preparing Data Object..")
  
  if(is.null(pThresh)){
    
    pThresh = 0.05
    
  }
  
  GMM <- vector(mode = "list", length = nrow(Table.ID))
  GMM.wFC <- vector(mode = "list", length = nrow(Table.ID))
  
  names(GMM) <- Table.ID[, 1]
  names(GMM.wFC) <- Table.ID[, 1]
  
  for(ii in 1:nrow(Table.ID)){
    
    GMM[[ii]] <- matrix(data = "NA", nrow = length(tableTopList), ncol = 4)
    colnames(GMM[[ii]]) <- c("Indiv", "clus", "FCvCaPval", "status")
    rownames(GMM[[ii]]) <- paste0("cond", 1:length(tableTopList))
    
    GMM.wFC[[ii]] <- matrix(data = "NA", nrow = length(tableTopList), ncol = 5)
    colnames(GMM.wFC[[ii]]) <- c("Indiv", "clus", "FCvCaPval", "status", "FC")
    rownames(GMM.wFC[[ii]]) <- paste0("cond", 1:length(tableTopList))
    
    for(jj in 1:length(tableTopList)){
      
      tt <- tableTopList[[jj]]
      
      if(names(GMM)[[ii]]%in%rownames(tt)){
        
        #
        GMM[[ii]][jj, 1] <- as.character(log(x = tt$adj.P.Val[which(rownames(tt)==names(GMM)[ii])]/pThresh, base = 2))
        GMM.wFC[[ii]][jj, 1] <- as.character(log(x = tt$adj.P.Val[which(rownames(tt)==names(GMM)[ii])]/pThresh, base = 2))
        
        if(!is.na(log(x = tt$adj.P.Val[which(rownames(tt)==names(GMM)[ii])]/pThresh, base = 2))){
          
          #
          if(tt$adj.P.Val[which(rownames(tt)==names(GMM)[ii])] < pThresh){
            
            GMM[[ii]][jj, 2] <- "P"
            GMM.wFC[[ii]][jj, 2] <- "P"
            
          } else {
            
            GMM[[ii]][jj, 2] <- "C"
            GMM.wFC[[ii]][jj, 2] <- "C"
            
          }
          
          #
          GMM[[ii]][jj, 3] <- as.character(tt$adj.P.Val[which(rownames(tt)==names(GMM)[ii])])
          GMM.wFC[[ii]][jj, 3] <- as.character(tt$adj.P.Val[which(rownames(tt)==names(GMM)[ii])])
          
          #
          GMM[[ii]][jj, 4] <- "OK"
          GMM.wFC[[ii]][jj, 4] <- "OK"
          
          #
          GMM.wFC[[ii]][jj, 5] <- as.character(tt$logFC[which(rownames(tt)==names(GMM)[ii])])
          
        }
        
      }
      
    }
    
  }
  
  dataInput <- list()
  dataInput[[1]] <- GMM.ID
  dataInput[[2]] <- GMM
  dataInput[[3]] <- GMM.wFC
  
  names(dataInput) <- c("IDmap", "res", "resFC")
  
  save(dataInput, file = "dataInput.RData")
  
  return(dataInput)
  
}