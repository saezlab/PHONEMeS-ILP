#
#  This file is part of the PHONEMeS-ILP method
#
#  Copyright (c) 2018 - RWTH Aachen - JRC COMBINE
#
#  File author(s): E.Gjerga (enio.gjerga@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  PHONEMeS website: https://saezlab.github.io/PHONEMeS/
#
##############################################################################
# $Id$

# Script showing how to build the input dataGMM object for PHONEMeS from the Limma
# output object

load(file = "allD_MOUSE.RData")
load(file = "ttop_list.RData")

tt1 <- ttop_list[[1]]
tt2 <- ttop_list[[2]]
tt3 <- ttop_list[[3]]
tt4 <- ttop_list[[4]]
tt5 <- ttop_list[[5]]
tt6 <- ttop_list[[6]]

###
measuredSpecies <- unique(c(tt1$ID, tt2$ID, tt3$ID, tt4$ID, tt5$ID, tt6$ID))
data.IDmap <- matrix(data = , nrow = length(measuredSpecies), ncol = 6)
colnames(data.IDmap) <- c("dataID", "UPID", "site", "res", "pos", "S.cc")
data.IDmap[, 1] <- measuredSpecies
data.IDmap[, 2] <- measuredSpecies
data.IDmap[, 6] <- measuredSpecies

data.IDmap <- as.data.frame(data.IDmap)
data.IDmap$dataID <- as.character(data.IDmap$dataID)
data.IDmap$UPID <- as.character(data.IDmap$UPID)
data.IDmap$S.cc <- as.character(data.IDmap$S.cc)

threshFC <- 0.05

##
v1.p <- rep(NA, nrow(data.IDmap))
v1.lo <- rep(NA, nrow(data.IDmap))
v1.fc <- rep(NA, nrow(data.IDmap))
v1.s <- rep(NA, nrow(data.IDmap))
v1.c <- rep(NA, nrow(data.IDmap))
for(i in 1:nrow(data.IDmap)){
  
  idx <- which(tt1$ID==data.IDmap$dataID[i])
  
  if(length(idx) > 0){
    
    v1.p[i] <- tt1$adj.P.Val[idx]
    v1.lo[i] <- log(x = v1.p[i]/threshFC)
    v1.fc[i] <- tt1$logFC[idx]
    v1.s[i] <- "OK"
    if(v1.lo[i] < 0){
      v1.c[i] <- "P"
    } else {
      v1.c[i] <- "C"
    }
    
  }
  
}

##
v2.p <- rep(NA, nrow(data.IDmap))
v2.lo <- rep(NA, nrow(data.IDmap))
v2.fc <- rep(NA, nrow(data.IDmap))
v2.s <- rep(NA, nrow(data.IDmap))
v2.c <- rep(NA, nrow(data.IDmap))
for(i in 1:nrow(data.IDmap)){
  
  idx <- which(tt2$ID==data.IDmap$dataID[i])
  
  if(length(idx) > 0){
    
    v2.p[i] <- tt2$adj.P.Val[idx]
    v2.lo[i] <- log(x = v2.p[i]/threshFC)
    v2.fc[i] <- tt2$logFC[idx]
    v2.s[i] <- "OK"
    if(v2.lo[i] < 0){
      v2.c[i] <- "P"
    } else {
      v2.c[i] <- "C"
    }
    
  }
  
}

##
v3.p <- rep(NA, nrow(data.IDmap))
v3.lo <- rep(NA, nrow(data.IDmap))
v3.fc <- rep(NA, nrow(data.IDmap))
v3.s <- rep(NA, nrow(data.IDmap))
v3.c <- rep(NA, nrow(data.IDmap))
for(i in 1:nrow(data.IDmap)){
  
  idx <- which(tt3$ID==data.IDmap$dataID[i])
  
  if(length(idx) > 0){
    
    v3.p[i] <- tt3$adj.P.Val[idx]
    v3.lo[i] <- log(x = v3.p[i]/threshFC)
    v3.fc[i] <- tt3$logFC[idx]
    v3.s[i] <- "OK"
    if(v3.lo[i] < 0){
      v3.c[i] <- "P"
    } else {
      v3.c[i] <- "C"
    }
    
  }
  
}

##
v4.p <- rep(NA, nrow(data.IDmap))
v4.lo <- rep(NA, nrow(data.IDmap))
v4.fc <- rep(NA, nrow(data.IDmap))
v4.s <- rep(NA, nrow(data.IDmap))
v4.c <- rep(NA, nrow(data.IDmap))
for(i in 1:nrow(data.IDmap)){
  
  idx <- which(tt4$ID==data.IDmap$dataID[i])
  
  if(length(idx) > 0){
    
    v4.p[i] <- tt4$adj.P.Val[idx]
    v4.lo[i] <- log(x = v4.p[i]/threshFC)
    v4.fc[i] <- tt4$logFC[idx]
    v4.s[i] <- "OK"
    if(v4.lo[i] < 0){
      v4.c[i] <- "P"
    } else {
      v4.c[i] <- "C"
    }
    
  }
  
}

##
v5.p <- rep(NA, nrow(data.IDmap))
v5.lo <- rep(NA, nrow(data.IDmap))
v5.fc <- rep(NA, nrow(data.IDmap))
v5.s <- rep(NA, nrow(data.IDmap))
v5.c <- rep(NA, nrow(data.IDmap))
for(i in 1:nrow(data.IDmap)){
  
  idx <- which(tt5$ID==data.IDmap$dataID[i])
  
  if(length(idx) > 0){
    
    v5.p[i] <- tt5$adj.P.Val[idx]
    v5.lo[i] <- log(x = v5.p[i]/threshFC)
    v5.fc[i] <- tt5$logFC[idx]
    v5.s[i] <- "OK"
    if(v5.lo[i] < 0){
      v5.c[i] <- "P"
    } else {
      v5.c[i] <- "C"
    }
    
  }
  
}

##
v6.p <- rep(NA, nrow(data.IDmap))
v6.lo <- rep(NA, nrow(data.IDmap))
v6.fc <- rep(NA, nrow(data.IDmap))
v6.s <- rep(NA, nrow(data.IDmap))
v6.c <- rep(NA, nrow(data.IDmap))
for(i in 1:nrow(data.IDmap)){
  
  idx <- which(tt6$ID==data.IDmap$dataID[i])
  
  if(length(idx) > 0){
    
    v6.p[i] <- tt6$adj.P.Val[idx]
    v6.lo[i] <- log(x = v6.p[i]/threshFC)
    v6.fc[i] <- tt6$logFC[idx]
    v6.s[i] <- "OK"
    if(v6.lo[i] < 0){
      v6.c[i] <- "P"
    } else {
      v6.c[i] <- "C"
    }
    
  }
  
}

##
GMM <- vector(mode = "list", length = nrow(data.IDmap))
GMM.wFC <- vector(mode = "list", length = nrow(data.IDmap))
names(GMM) <- data.IDmap$dataID
names(GMM.wFC) <- data.IDmap$dataID
for(i in 1:length(GMM)){
  
  GMM[[i]]<-rbind(
    c(as.character(v1.lo[i]), as.character(v1.c[i]), as.character(v1.p[i]), as.character(v1.s[i])),
    c(as.character(v2.lo[i]), as.character(v2.c[i]), as.character(v2.p[i]), as.character(v2.s[i])),
    c(as.character(v3.lo[i]), as.character(v3.c[i]), as.character(v3.p[i]), as.character(v3.s[i])),
    c(as.character(v4.lo[i]), as.character(v4.c[i]), as.character(v4.p[i]), as.character(v4.s[i])),
    c(as.character(v5.lo[i]), as.character(v5.c[i]), as.character(v5.p[i]), as.character(v5.s[i])),
    c(as.character(v6.lo[i]), as.character(v6.c[i]), as.character(v6.p[i]), as.character(v6.s[i]))
  )
  colnames(GMM[[i]])<-c("Indiv", "clus","FCvCaPval","status")
  rownames(GMM[[i]])<-names(ttop_list)
  
  GMM.wFC[[i]]<-rbind(
    c(as.character(v1.lo[i]), as.character(v1.c[i]), as.character(v1.p[i]), as.character(v1.s[i]), as.character(v1.fc[i])),
    c(as.character(v2.lo[i]), as.character(v2.c[i]), as.character(v2.p[i]), as.character(v2.s[i]), as.character(v2.fc[i])),
    c(as.character(v3.lo[i]), as.character(v3.c[i]), as.character(v3.p[i]), as.character(v3.s[i]), as.character(v3.fc[i])),
    c(as.character(v4.lo[i]), as.character(v4.c[i]), as.character(v4.p[i]), as.character(v4.s[i]), as.character(v4.fc[i])),
    c(as.character(v5.lo[i]), as.character(v5.c[i]), as.character(v5.p[i]), as.character(v5.s[i]), as.character(v5.fc[i])),
    c(as.character(v6.lo[i]), as.character(v6.c[i]), as.character(v6.p[i]), as.character(v6.s[i]), as.character(v6.fc[i]))
  )
  colnames(GMM.wFC[[i]])<-c("Indiv", "clus", "FCvCaPval", "status", "logFC")
  rownames(GMM.wFC[[i]])<-names(ttop_list)
  
}




###
# Saving GMM object as a list
data.IDmap <- data.IDmap[, c(1, 2, 6)]
GMM.ID <- data.IDmap
colnames(GMM.ID) <- c("dataID", "UPID", "S.cc")
GMM.ID <- as.data.frame(GMM.ID)

toGroup <- c("JAK", "STAT", "PIK3C", "AKT", "PTPN", "GAB", "ERK")

#
meas <- names(GMM)[which(grepl(pattern = "JAK", x = names(GMM)))]
if(length(meas) > 0){
  idx <- which(names(GMM)%in%meas)
  for(i in 1:length(idx)){
    names(GMM)[idx[i]] <- paste0("JAK_", strsplit(x = names(GMM)[idx[i]], split = "_")[[1]][[2]])
  }
}

#
meas <- names(GMM)[which(grepl(pattern = "STAT", x = names(GMM)))]
if(length(meas) > 0){
  idx <- which(names(GMM)%in%meas)
  for(i in 1:length(idx)){
    names(GMM)[idx[i]] <- paste0("STAT_", strsplit(x = names(GMM)[idx[i]], split = "_")[[1]][[2]])
  }
}

#
meas <- names(GMM)[which(grepl(pattern = "PIK3", x = names(GMM)))]
if(length(meas) > 0){
  idx <- which(names(GMM)%in%meas)
  for(i in 1:length(idx)){
    names(GMM)[idx[i]] <- paste0("PIK3C_", strsplit(x = names(GMM)[idx[i]], split = "_")[[1]][[2]])
  }
}

#
meas <- names(GMM)[which(grepl(pattern = "AKT", x = names(GMM)))]
meas <- meas[-which(grepl(pattern = "AKT1S", x = meas))]
meas <- meas[-which(grepl(pattern = "AKTIP", x = meas))]
if(length(meas) > 0){
  idx <- which(names(GMM)%in%meas)
  for(i in 1:length(idx)){
    names(GMM)[idx[i]] <- paste0("AKT_", strsplit(x = names(GMM)[idx[i]], split = "_")[[1]][[2]])
  }
}

#
meas <- names(GMM)[which(grepl(pattern = "PTPN", x = names(GMM)))]
if(length(meas) > 0){
  idx <- which(names(GMM)%in%meas)
  for(i in 1:length(idx)){
    names(GMM)[idx[i]] <- paste0("PTPN_", strsplit(x = names(GMM)[idx[i]], split = "_")[[1]][[2]])
  }
}

#
meas <- names(GMM)[which(grepl(pattern = "GAB", x = names(GMM)))]
meas <- meas[-which(grepl(pattern = "GABP", x = meas))]
if(length(meas) > 0){
  idx <- which(names(GMM)%in%meas)
  for(i in 1:length(idx)){
    names(GMM)[idx[i]] <- paste0("GAB_", strsplit(x = names(GMM)[idx[i]], split = "_")[[1]][[2]])
  }
}

#
meas <- c(names(GMM)[which(grepl(pattern = "MAPK1_", x = names(GMM)))], names(GMM)[which(grepl(pattern = "MAPK3_", x = names(GMM)))])
if(length(meas) > 0){
  idx <- which(names(GMM)%in%meas)
  for(i in 1:length(idx)){
    names(GMM)[idx[i]] <- paste0("ERK_", strsplit(x = names(GMM)[idx[i]], split = "_")[[1]][[2]])
  }
}

names(GMM.wFC) <- names(GMM)

###
# GMM.ID
meas <- GMM.ID$dataID[which(grepl(pattern = "JAK", x = GMM.ID$dataID))]
if(length(meas) > 0){
  idx <- which(GMM.ID$dataID%in%meas)
  for(i in 1:length(idx)){
    GMM.ID$dataID[idx[i]] <- paste0("JAK_", strsplit(x = GMM.ID$dataID[idx[i]], split = "_")[[1]][[2]])
  }
}

meas <- GMM.ID$dataID[which(grepl(pattern = "STAT", x = GMM.ID$dataID))]
if(length(meas) > 0){
  idx <- which(GMM.ID$dataID%in%meas)
  for(i in 1:length(idx)){
    GMM.ID$dataID[idx[i]] <- paste0("STAT_", strsplit(x = GMM.ID$dataID[idx[i]], split = "_")[[1]][[2]])
  }
}

meas <- GMM.ID$dataID[which(grepl(pattern = "AKT", x = GMM.ID$dataID))]
meas <- meas[-which(grepl(pattern = "AKT1S", x = meas))]
meas <- meas[-which(grepl(pattern = "AKTIP", x = meas))]
if(length(meas) > 0){
  idx <- which(GMM.ID$dataID%in%meas)
  for(i in 1:length(idx)){
    GMM.ID$dataID[idx[i]] <- paste0("AKT_", strsplit(x = GMM.ID$dataID[idx[i]], split = "_")[[1]][[2]])
  }
}

meas <- GMM.ID$dataID[which(grepl(pattern = "PTPN", x = GMM.ID$dataID))]
if(length(meas) > 0){
  idx <- which(GMM.ID$dataID%in%meas)
  for(i in 1:length(idx)){
    GMM.ID$dataID[idx[i]] <- paste0("PTPN_", strsplit(x = GMM.ID$dataID[idx[i]], split = "_")[[1]][[2]])
  }
}

meas <- GMM.ID$dataID[which(grepl(pattern = "GAB", x = GMM.ID$dataID))]
meas <- meas[-which(grepl(pattern = "GABP", x = meas))]
if(length(meas) > 0){
  idx <- which(GMM.ID$dataID%in%meas)
  for(i in 1:length(idx)){
    GMM.ID$dataID[idx[i]] <- paste0("GAB_", strsplit(x = GMM.ID$dataID[idx[i]], split = "_")[[1]][[2]])
  }
}

meas <- c(GMM.ID$dataID[which(grepl(pattern = "MAPK1_", x = GMM.ID$dataID))], GMM.ID$dataID[which(grepl(pattern = "MAPK3_", x = GMM.ID$dataID))])
if(length(meas) > 0){
  idx <- which(GMM.ID$dataID%in%meas)
  for(i in 1:length(idx)){
    GMM.ID$dataID[idx[i]] <- paste0("ERK_", strsplit(x = GMM.ID$dataID[idx[i]], split = "_")[[1]][[2]])
  }
}

GMM.ID$UPID <- GMM.ID$dataID
GMM.ID$S.cc <- GMM.ID$dataID


save(list=c("GMM.ID", "GMM","GMM.wFC"), file="dataGMM.RData")
