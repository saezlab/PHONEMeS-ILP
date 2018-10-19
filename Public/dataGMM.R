library(readr)
library(dplyr)
library(data.table)
library(tidyr)

#' Generate GMMres object
#' 
#' @param df DataFrame containing phosphosite measurement information.
#' @param threshold When to consider a phosphosite as perturbed.
#
#' @return GMMres class object.
getGMMres <- function(df, threshold=0.1){
  condition_list <- df %>% pull(term) %>% unique()
  
  condition_data1 <- df %>% filter(term==condition_list[1])
  
  # define mapping table
  mappingTable <- matrix(data = , nrow = nrow(condition_data1), ncol = 6)
  colnames(mappingTable) <- c("dataID", "UPID", "site", "res", "pos", "S.cc")
  mappingTable[, 1] <- condition_data1$gene
  for(i in 1:nrow(condition_data1)){
    mappingTable[i, 2] <- strsplit(x = condition_data1$ProteinName[i], split = "|", fixed = TRUE)[[1]][length(strsplit(x = condition_data1$ProteinName[i], split = "|", fixed = TRUE)[[1]])]
    mappingTable[i, 3] <- condition_data1$residues_str[i]
    mappingTable[i, 4] <- substr(x = condition_data1$residues_str[i], start = 1, stop = 1)
    mappingTable[i, 5] <- condition_data1$pep_pos[i]
    mappingTable[i, 6] <- paste0(mappingTable[i, 2], "_", mappingTable[i, 3])
  }
  
  GMM.ID <- as.data.frame(mappingTable)
  GMM.ID$dataID <- as.character(GMM.ID$dataID)
  GMM.ID$UPID <- as.character(GMM.ID$UPID)
  GMM.ID$site <- as.character(GMM.ID$site)
  GMM.ID$pos <- as.character(GMM.ID$pos)
  GMM.ID$S.cc <- as.character(GMM.ID$S.cc)
  
  # prepare data
  v.p <- list()
  v.lo <- list()
  v.fc <- list()
  v.s <- list()
  v.c <- list()
  for (i_condition in condition_list){
    i_condition_data <- df %>% filter(term==i_condition)
    v.p[[i_condition]] <- rep(NA, nrow(GMM.ID))
    v.lo[[i_condition]] <- rep(NA, nrow(GMM.ID))
    v.fc[[i_condition]] <- rep(NA, nrow(GMM.ID))
    v.s[[i_condition]] <- rep(NA, nrow(GMM.ID))
    v.c[[i_condition]] <- rep(NA, nrow(GMM.ID))
    for(i in 1:nrow(GMM.ID)){
      v.p[[i_condition]][i] <- i_condition_data$p.value[i]
      v.lo[[i_condition]][i] <- log2(x = v.p[[i_condition]][i]/threshold)
      v.fc[[i_condition]][i] <- i_condition_data$estimate[i]
      v.s[[i_condition]][i] <- "OK"
      if(v.lo[[i_condition]][i]<0){
        v.c[[i_condition]][i] <- "P"
      } else {
        v.c[[i_condition]][i] <- "C"
      }
    }
  }
  
  # Prepare GMM
  GMM <- vector(mode = "list", length = nrow(GMM.ID))
  GMM.wFC <- vector(mode = "list", length = nrow(GMM.ID))
  names(GMM) <- GMM.ID$dataID
  names(GMM.wFC) <- GMM.ID$dataID
  for(i in 1:length(GMM)){
    # initialize with first condition
    GMM[[i]] <- c(as.character(v.lo[[condition_list[1]]][i]), 
                  as.character(v.c[[condition_list[1]]][i]), 
                  as.character(v.p[[condition_list[1]]][i]), 
                  as.character(v.s[[condition_list[1]]][i]))
    GMM.wFC[[i]] <- c(as.character(v.lo[[condition_list[1]]][i]), 
                      as.character(v.c[[condition_list[1]]][i]), 
                      as.character(v.p[[condition_list[1]]][i]), 
                      as.character(v.s[[condition_list[1]]][i]), 
                      as.character(v.fc[[condition_list[1]]][i]))
    
    # include information from the rest of conditions (if any)
    for (j in 2:length(condition_list)){
      GMM_to_add <- c(as.character(v.lo[[condition_list[j]]][i]), 
                      as.character(v.c[[condition_list[j]]][i]), 
                      as.character(v.p[[condition_list[j]]][i]), 
                      as.character(v.s[[condition_list[j]]][i]))
      GMM[[i]]<-rbind(GMM[[i]], GMM_to_add)
      
      GMM.wFC_to_add <- c(as.character(v.lo[[condition_list[j]]][i]), 
                          as.character(v.c[[condition_list[j]]][i]), 
                          as.character(v.p[[condition_list[j]]][i]), 
                          as.character(v.s[[condition_list[j]]][i]), 
                          as.character(v.fc[[condition_list[j]]][i]))
      GMM.wFC[[i]]<-rbind(GMM.wFC[[i]], GMM.wFC_to_add)
    }
    
    # set column and row names
    colnames(GMM[[i]]) <- c("Indiv", "clus","FCvCaPval","status")
    rownames(GMM[[i]]) <- condition_list
    colnames(GMM.wFC[[i]]) <- c("Indiv", "clus", "FCvCaPval", "status", "logFC")
    rownames(GMM.wFC[[i]]) <- condition_list
  }
  
  ###
  # GMM object as a list
  GMM.ID <- GMM.ID[, c(1, 2, 6)]
  GMM.ID <- GMM.ID
  colnames(GMM.ID) <- c("dataID", "UPID", "S.cc")
  GMM.ID <- as.data.frame(GMM.ID)
  
  # GMMres class object
  dataGMM<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)
  
  return(dataGMM)
}
