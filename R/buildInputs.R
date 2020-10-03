#' Build Inputs for PHONEMeS
#' 
#' @param tableTopList 
#' @param fcThresh
#' @param pThresh
#' @param idxID
#' @param idxFC
#' @param idxPval
#' @param mappingTable
#' @param namesConditions
#' @param directionality
#' @param solver Solver to use for solving the ILP.
#
#' @return A data-input object for PHONEMeS.

buildInputs <- function(tableTopList = tableTopList, 
                        fcThresh = NULL, 
                        pThresh = 0.05, 
                        idxID = 1, 
                        idxFC = 2, 
                        idxPval = 6, 
                        mappingTable = NULL, 
                        namesConditions = NULL, 
                        directionality = "no"){
  
  ###
  # Build GMM.ID object
  allowed_directionality <- c("up", "down", "no")
  if((directionality%in%allowed_directionality)==FALSE){
    stop("Wrong directionality value. It should either be no/up/down")
  }
  
  if(is.null(mappingTable)){
    
    species <- c()
    for(i in 1:length(tableTopList)){
      
      species <- c(species, tableTopList[[i]][, idxID])
      
    }
    
    species <- unique(species)
    
    proteins <- unlist(lapply(strsplit(x = species, split = "_", fixed = TRUE), '[[', 1))
    residues <- unlist(lapply(strsplit(x = species, split = "_", fixed = TRUE), '[[', 2))
    
    if(length(which(residues==""))){
      proteins <- proteins[-which(residues=="")]
      residues <- residues[-which(residues=="")]
    }
    
    Table.ID <- matrix(data = , nrow = length(proteins), ncol = 6)
    colnames(Table.ID) <- c("dataID", "UPID", "site", "pos", "res", "S.cc")
    
    Table.ID[, 1] <- paste0(proteins, "_", residues)
    Table.ID[, 3] <- residues
    for(i in 1:nrow(Table.ID)){
      
      Table.ID[i, 2] <- strsplit(x = Table.ID[i, 1], split = "_", fixed = TRUE)[[1]][1]
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
    
    pThresh = rep(0.05, length(ttop_list))
    
  }
  
  GMM <- vector(mode = "list", length = nrow(Table.ID))
  GMM.wFC <- vector(mode = "list", length = nrow(Table.ID))
  
  names(GMM) <- Table.ID[, 1]
  names(GMM.wFC) <- Table.ID[, 1]
  
  for(ii in 1:nrow(Table.ID)){
    
    GMM[[ii]] <- matrix(data = "NA", nrow = length(tableTopList), ncol = 4)
    colnames(GMM[[ii]]) <- c("Indiv", "clus", "FCvCaPval", "status")
    if((is.null(namesConditions)) || (length(namesConditions)!=length(tableTopList))){
      rownames(GMM[[ii]]) <- paste0("cond", 1:length(tableTopList))
    } else {
      rownames(GMM[[ii]]) <- namesConditions
    }
    
    GMM.wFC[[ii]] <- matrix(data = "NA", nrow = length(tableTopList), ncol = 5)
    colnames(GMM.wFC[[ii]]) <- c("Indiv", "clus", "FCvCaPval", "status", "FC")
    if((is.null(namesConditions)) || (length(namesConditions)!=length(tableTopList))){
      rownames(GMM.wFC[[ii]]) <- paste0("cond", 1:length(tableTopList))
    } else {
      rownames(GMM.wFC[[ii]]) <- namesConditions
    }
    
    for(jj in 1:length(tableTopList)){
      
      tt <- tableTopList[[jj]]
      
      if(names(GMM)[[ii]]%in%rownames(tt)){
        
        #
        GMM[[ii]][jj, 1] <- as.character(log(x = tt[which(rownames(tt)==names(GMM)[ii]), idxPval]/pThresh[jj], base = 2))
        GMM.wFC[[ii]][jj, 1] <- as.character(log(x = tt[which(rownames(tt)==names(GMM)[ii]), idxPval]/pThresh[jj], base = 2))
        
        if(!is.na(log(x = tt[which(rownames(tt)==names(GMM)[ii]), idxPval]/pThresh[jj], base = 2))){
          
          #
          if(is.null(fcThresh)){
            if(tt[which(rownames(tt)==names(GMM)[ii]), idxPval] <= pThresh[jj]){
              
              GMM[[ii]][jj, 2] <- "P"
              GMM.wFC[[ii]][jj, 2] <- "P"
              
            } else {
              
              GMM[[ii]][jj, 2] <- "C"
              GMM.wFC[[ii]][jj, 2] <- "C"
              
            }
          } else {
            
            if(directionality=="no"){
              
              if(abs(tt[which(rownames(tt)==names(GMM)[ii]), idxFC]) >= fcThresh[jj]){
                
                GMM[[ii]][jj, 2] <- "P"
                GMM.wFC[[ii]][jj, 2] <- "P"
                
              } else {
                
                GMM[[ii]][jj, 2] <- "C"
                GMM.wFC[[ii]][jj, 2] <- "C"
                
              }
              
            } else {
              
              if(directionality=="up"){
                
                if(tt[which(rownames(tt)==names(GMM)[ii]), idxFC] >= fcThresh[jj]){
                  
                  GMM[[ii]][jj, 2] <- "P"
                  GMM.wFC[[ii]][jj, 2] <- "P"
                  
                } else {
                  
                  GMM[[ii]][jj, 2] <- "C"
                  GMM.wFC[[ii]][jj, 2] <- "C"
                  
                }
                
              } else {
                
                if(tt[which(rownames(tt)==names(GMM)[ii]), idxFC] <= fcThresh[jj]){
                  
                  GMM[[ii]][jj, 2] <- "P"
                  GMM.wFC[[ii]][jj, 2] <- "P"
                  
                } else {
                  
                  GMM[[ii]][jj, 2] <- "C"
                  GMM.wFC[[ii]][jj, 2] <- "C"
                  
                }
                
              }
            }
          }
          
          #
          GMM[[ii]][jj, 3] <- as.character(tt[which(rownames(tt)==names(GMM)[ii]), idxPval])
          GMM.wFC[[ii]][jj, 3] <- as.character(tt[which(rownames(tt)==names(GMM)[ii]), idxPval])
          
          #
          GMM[[ii]][jj, 4] <- "OK"
          GMM.wFC[[ii]][jj, 4] <- "OK"
          
          #
          GMM.wFC[[ii]][jj, 5] <- as.character(tt[which(rownames(tt)==names(GMM)[ii]), idxFC])
          
        }
        
      }
      
    }
    
  }
  
  for(ii in 1:length(GMM)){
    
    for(jj in 1:nrow(GMM[[ii]])){
      
      if((GMM[[ii]][jj, 2]=="C") && (as.numeric(GMM[[ii]][jj, 1])<0)){
        
        GMM[[ii]][jj, 1] <- as.character((-1)*(as.numeric(GMM[[ii]][jj, 1])))
        GMM.wFC[[ii]][jj, 1] <- as.character((-1)*(as.numeric(GMM.wFC[[ii]][jj, 1])))
                                             
      }
      
    }
    
  }
  
  dataInput <- list()
  dataInput[[1]] <- GMM.ID
  dataInput[[2]] <- GMM
  dataInput[[3]] <- GMM.wFC
  
  names(dataInput) <- c("IDmap", "res", "resFC")
  
  dataInput<-new("GMMres", res=dataInput$res, IDmap=dataInput$IDmap, resFC=dataInput$resFC)
  
  return(dataInput)
  
}
