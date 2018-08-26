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

# Script showing how to build the background network object to train in PHONEMeS

library(readr)

omnipath_mouse <- read_delim("omnipath_mouse.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
omnipath_mouse_enzyme_substrate_all <- read_delim("omnipath_mouse_enzyme-substrate_all.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#
idx1 <- intersect(x = which(!is.na(omnipath_mouse$`Direction_A-B`)), y = which(!is.na(omnipath_mouse$`Stimulatory_A-B`)))
idx2 <- intersect(x = which(!is.na(omnipath_mouse$`Direction_A-B`)), y = which(!is.na(omnipath_mouse$`Inhibitory_A-B`)))
idx <- unique(c(idx1, idx2))

ppi <- omnipath_mouse[idx, c(1, 2, 3, 4)]
ppi$GeneSymbol_A <- toupper(x = ppi$GeneSymbol_A)
ppi$GeneSymbol_B <- toupper(x = ppi$GeneSymbol_B)

#
mapping <- matrix(data = , nrow = 2*nrow(ppi), ncol = 2)
mapping[1:nrow(ppi), 1] <- ppi$UniProt_A
mapping[(nrow(ppi)+1):(2*nrow(ppi)), 1] <- ppi$UniProt_B
mapping[1:nrow(ppi), 2] <- ppi$GeneSymbol_A
mapping[(nrow(ppi)+1):(2*nrow(ppi)), 2] <- ppi$GeneSymbol_B
mapping <- unique(mapping)
colnames(mapping) <- c("Uniprot", "Gene")

ppi <- ppi[, c(2, 4)]
idx1 <- which(grepl(pattern = "PDGFR", x = ppi$GeneSymbol_A))
idx2 <- which(grepl(pattern = "PDGFR", x = ppi$GeneSymbol_B))
idx <- unique(c(idx1, idx2))
ppi <- ppi[-idx, ]

idx1 <- which(omnipath_mouse_enzyme_substrate_all$modification=="phosphorylation")
idx2 <- which(omnipath_mouse_enzyme_substrate_all$modification=="dephosphorylation")

idx <- c(idx1, idx2)

omnipath_mouse_enzyme_substrate_all <- omnipath_mouse_enzyme_substrate_all[idx, ]
ptms <- matrix(data = , nrow = nrow(omnipath_mouse_enzyme_substrate_all), ncol = 2)
for(ii in 1:nrow(omnipath_mouse_enzyme_substrate_all)){
  
  idx1 <- which(mapping[, 1]==omnipath_mouse_enzyme_substrate_all$enzyme[ii])
  idx2 <- which(mapping[, 1]==omnipath_mouse_enzyme_substrate_all$substrate[ii])
  
  if(length(idx1)>0){
    ptms[ii, 1] <- mapping[idx1[1], 2]
  } else {
    ptms[ii, 1] <- omnipath_mouse_enzyme_substrate_all$enzyme[ii]
  }
  
  if(length(idx2)>0){
    ptms[ii, 2] <- paste0(mapping[idx2[1], 2], "_", omnipath_mouse_enzyme_substrate_all$residue[ii], omnipath_mouse_enzyme_substrate_all$offset[ii])
  } else {
    ptms[ii, 2] <- paste0(omnipath_mouse_enzyme_substrate_all$substrate[ii], "_", omnipath_mouse_enzyme_substrate_all$residue[ii], omnipath_mouse_enzyme_substrate_all$offset[ii])
  }
  
}
colnames(ptms) <- c("Kinase", "Substrate")

#
string_interactions_pdgfra <- read_delim("string_interactions_pdgfra.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
string_interactions_pdgfra <- string_interactions_pdgfra[1:nrow(string_interactions_pdgfra), c(1, 2)]
string_interactions_pdgfra[, 1] <- toupper(as.matrix(string_interactions_pdgfra[, 1]))
string_interactions_pdgfra[, 2] <- toupper(as.matrix(string_interactions_pdgfra[, 2]))
string_interactions_pdgfra <- string_interactions_pdgfra[which(string_interactions_pdgfra[, 1]=="PDGFA"), ]
string_interactions_pdgfra <- string_interactions_pdgfra[4:8, ]
colnames(string_interactions_pdgfra) <- colnames(ppi)
ppi <- rbind(ppi, string_interactions_pdgfra)

#
string_interactions_pdgfra <- read_delim("string_interactions_pdgfrb.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
string_interactions_pdgfra <- string_interactions_pdgfra[1:nrow(string_interactions_pdgfra), c(1, 2)]
string_interactions_pdgfra[, 1] <- toupper(as.matrix(string_interactions_pdgfra[, 1]))
string_interactions_pdgfra[, 2] <- toupper(as.matrix(string_interactions_pdgfra[, 2]))
string_interactions_pdgfra <- string_interactions_pdgfra[which(grepl(pattern = "PDGFR", x = as.matrix(string_interactions_pdgfra[, 1]))), ]
string_interactions_pdgfra <- string_interactions_pdgfra[c(4,8,9,10), ]
colnames(string_interactions_pdgfra) <- colnames(ppi)
ppi <- rbind(ppi, string_interactions_pdgfra)

ppi <- unique(ppi)

geneset_pdgf <- read_csv("geneset_pdgf.txt", col_names = FALSE)
genes <- geneset_pdgf$X1
genes[which(grepl(pattern = "PDGF", x = genes))]

idx1 <- which(ppi$GeneSymbol_A%in%genes)
idx2 <- which(ppi$GeneSymbol_B%in%genes)

idx <- intersect(x = idx1, y = idx2)

ppi <- ppi[idx, ]

###
BN <- matrix(data = , nrow = nrow(ptms), ncol = 8)
colnames(BN) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

cnt <- 1
for(i in 1:nrow(ptms)){
  
  BN[i, 1] <- strsplit(x = ptms[i, 2], split = "_")[[1]][[1]]
  BN[i, 2] <- strsplit(x = ptms[i, 2], split = "_")[[1]][[1]]
  BN[i, 3] <- ptms[i, 1]
  BN[i, 4] <- ptms[i, 1]
  BN[i, 5] <- substr(x = strsplit(x = ptms[i, 2], split = "_")[[1]][[2]], start = 1, stop = 1)
  BN[i, 6] <- substr(x = strsplit(x = ptms[i, 2], split = "_")[[1]][[2]], start = 2, stop = nchar(strsplit(x = ptms[i, 2], split = "_")[[1]][[2]]))
  BN[i, 7] <- paste0("e", cnt)
  cnt <- cnt+1
  BN[i, 8] <- ptms[i, 2]
  
}

BN1 <- BN

###
BN <- matrix(data = , nrow = nrow(ppi), ncol = 8)
colnames(BN) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

for(i in 1:nrow(ppi)){
  
  BN[i, 1] <- as.character(ppi[i, 2])
  BN[i, 2] <- as.character(ppi[i, 2])
  BN[i, 3] <- as.character(ppi[i, 1])
  BN[i, 4] <- as.character(ppi[i, 1])
  BN[i, 5] <- "R"
  BN[i, 6] <- "1"
  BN[i, 7] <- paste0("e", cnt)
  cnt <- cnt+1
  BN[i, 8] <- paste0(BN[i, 1], "_R1")
  
}

BN2 <- BN

###
BN <- rbind(BN1, BN2)

BN <- as.data.frame(BN)
BN$S.AC <- as.character(BN$S.AC)
BN$S.ID <- as.character(BN$S.ID)
BN$K.AC <- as.character(BN$K.AC)
BN$K.ID <- as.character(BN$K.ID)
BN$res <- as.character(BN$res)
BN$pos <- as.character(BN$pos)
BN$SID <- as.character(BN$SID)
BN$S.cc <- as.character(BN$S.cc)

##
toGroup <- c("JAK", "STAT", "PIK3C", "AKT", "PTPN", "GAB", "ERK")

#
allJAK <- unique(c(BN$S.ID[which(grepl(pattern = "JAK", x = BN$S.ID))], BN$K.ID[which(grepl(pattern = "JAK", x = BN$K.ID))], BN$S.cc[which(grepl(pattern = "JAK", x = BN$S.cc))]))

idx1 <- which(BN$S.ID%in%allJAK)
if(length(idx1) > 0){BN$S.ID[idx1] <- "JAK"}

idx2 <- which(BN$K.ID%in%allJAK)
if(length(idx2)>0){BN$K.ID[idx2] <- "JAK"}

idx3 <- which(BN$S.cc%in%allJAK)
if(length(idx3)>0){for(i in 1:length(idx3)){BN$S.cc[idx3[i]]<-paste0("JAK_", strsplit(x = BN$S.cc[idx3[i]], split = "_")[[1]][2])}}

uniprot <- BN$S.AC[which(BN$S.ID=="JAK")[1]]

BN$S.AC[idx1] <- "JAK"

BN$K.AC[idx2] <- "JAK"

#
allSTAT <- unique(c(BN$S.ID[which(grepl(pattern = "STAT", x = BN$S.ID))], BN$K.ID[which(grepl(pattern = "STAT", x = BN$K.ID))], BN$S.cc[which(grepl(pattern = "STAT", x = BN$S.cc))]))

idx1 <- which(BN$S.ID%in%allSTAT)
if(length(idx1) > 0){BN$S.ID[idx1] <- "STAT"}

idx2 <- which(BN$K.ID%in%allSTAT)
if(length(idx2)>0){BN$K.ID[idx2] <- "STAT"}

idx3 <- which(BN$S.cc%in%allSTAT)
if(length(idx3)>0){for(i in 1:length(idx3)){BN$S.cc[idx3[i]]<-paste0("STAT_", strsplit(x = BN$S.cc[idx3[i]], split = "_")[[1]][2])}}

uniprot <- BN$S.AC[which(BN$S.ID=="STAT")[1]]

BN$S.AC[idx1] <- "STAT"

BN$K.AC[idx2] <- "STAT"

#
allPIK3C <- unique(c(BN$S.ID[which(grepl(pattern = "PIK3", x = BN$S.ID))], BN$K.ID[which(grepl(pattern = "PIK3", x = BN$K.ID))], BN$S.cc[which(grepl(pattern = "PIK3", x = BN$S.cc))]))

idx1 <- which(BN$S.ID%in%allPIK3C)
if(length(idx1) > 0){BN$S.ID[idx1] <- "PIK3C"}

idx2 <- which(BN$K.ID%in%allPIK3C)
if(length(idx2)>0){BN$K.ID[idx2] <- "PIK3C"}

idx3 <- which(BN$S.cc%in%allPIK3C)
if(length(idx3)>0){for(i in 1:length(idx3)){BN$S.cc[idx3[i]]<-paste0("PIK3C_", strsplit(x = BN$S.cc[idx3[i]], split = "_")[[1]][2])}}

uniprot <- BN$S.AC[which(BN$S.ID=="PIK3C")[1]]

BN$S.AC[idx1] <- "PIK3C"

BN$K.AC[idx2] <- "PIK3C"

#
allAKT <- unique(c(BN$S.ID[which(grepl(pattern = "AKT", x = BN$S.ID))], BN$K.ID[which(grepl(pattern = "AKT", x = BN$K.ID))], BN$S.cc[which(grepl(pattern = "AKT", x = BN$S.cc))]))
allAKT <- allAKT[-which(grepl(pattern = "AKT1S1", x = allAKT))]

idx1 <- which(BN$S.ID%in%allAKT)
if(length(idx1) > 0){BN$S.ID[idx1] <- "AKT"}

idx2 <- which(BN$K.ID%in%allAKT)
if(length(idx2)>0){BN$K.ID[idx2] <- "AKT"}

idx3 <- which(BN$S.cc%in%allAKT)
if(length(idx3)>0){for(i in 1:length(idx3)){BN$S.cc[idx3[i]]<-paste0("AKT_", strsplit(x = BN$S.cc[idx3[i]], split = "_")[[1]][2])}}

uniprot <- BN$S.AC[which(BN$S.ID=="AKT")[1]]

BN$S.AC[idx1] <- "AKT"

BN$K.AC[idx2] <- "AKT"

#
allPTPN <- unique(c(BN$S.ID[which(grepl(pattern = "PTPN", x = BN$S.ID))], BN$K.ID[which(grepl(pattern = "PTPN", x = BN$K.ID))], BN$S.cc[which(grepl(pattern = "PTPN", x = BN$S.cc))]))

idx1 <- which(BN$S.ID%in%allPTPN)
if(length(idx1) > 0){BN$S.ID[idx1] <- "PTPN"}

idx2 <- which(BN$K.ID%in%allPTPN)
if(length(idx2)>0){BN$K.ID[idx2] <- "PTPN"}

idx3 <- which(BN$S.cc%in%allPTPN)
if(length(idx3)>0){for(i in 1:length(idx3)){BN$S.cc[idx3[i]]<-paste0("PTPN_", strsplit(x = BN$S.cc[idx3[i]], split = "_")[[1]][2])}}

uniprot <- BN$S.AC[which(BN$S.ID=="PTPN")[1]]

BN$S.AC[idx1] <- "PTPN"

BN$K.AC[idx2] <- "PTPN"

#
allGAB <- unique(c(BN$S.ID[which(grepl(pattern = "GAB", x = BN$S.ID))], BN$K.ID[which(grepl(pattern = "GAB", x = BN$K.ID))], BN$S.cc[which(grepl(pattern = "GAB", x = BN$S.cc))]))
allGAB <- allGAB[-which(grepl(pattern = "GABB", x = allGAB))]

idx1 <- which(BN$S.ID%in%allGAB)
if(length(idx1) > 0){BN$S.ID[idx1] <- "GAB"}

idx2 <- which(BN$K.ID%in%allGAB)
if(length(idx2)>0){BN$K.ID[idx2] <- "GAB"}

idx3 <- which(BN$S.cc%in%allGAB)
if(length(idx3)>0){for(i in 1:length(idx3)){BN$S.cc[idx3[i]]<-paste0("GAB_", strsplit(x = BN$S.cc[idx3[i]], split = "_")[[1]][2])}}

uniprot <- BN$S.AC[which(BN$S.ID=="GAB")[1]]

BN$S.AC[idx1] <- "GAB"

BN$K.AC[idx2] <- "GAB"

#
allERK <- unique(c(BN$S.ID[which(grepl(pattern = "MAPK", x = BN$S.ID))], BN$K.ID[which(grepl(pattern = "MAPK", x = BN$K.ID))], BN$S.cc[which(grepl(pattern = "MAPK", x = BN$S.cc))]))
allERK <- allERK[c(which(grepl(pattern = "MAPK1", x = allERK)), which(grepl(pattern = "MAPK3", x = allERK)))]
allERK <- allERK[-c(which(grepl(pattern = "MAPK10", x = allERK)),which(grepl(pattern = "MAPK11", x = allERK)),which(grepl(pattern = "MAPK12", x = allERK)),
                    which(grepl(pattern = "MAPK13", x = allERK)),which(grepl(pattern = "MAPK15", x = allERK)),which(grepl(pattern = "MAPK14", x = allERK)))]

idx1 <- which(BN$S.ID%in%allERK)
if(length(idx1) > 0){BN$S.ID[idx1] <- "ERK"}

idx2 <- which(BN$K.ID%in%allERK)
if(length(idx2)>0){BN$K.ID[idx2] <- "ERK"}

idx3 <- which(BN$S.cc%in%allERK)
if(length(idx3)>0){for(i in 1:length(idx3)){BN$S.cc[idx3[i]]<-paste0("ERK_", strsplit(x = BN$S.cc[idx3[i]], split = "_")[[1]][2])}}

uniprot <- BN$S.AC[which(BN$S.ID=="ERK")[1]]

BN$S.AC[idx1] <- "ERK"

BN$K.AC[idx2] <- "ERK"

#
allPDGF <- unique(c(BN$S.ID[which(grepl(pattern = "PDGF", x = BN$S.ID))], BN$K.ID[which(grepl(pattern = "PDGF", x = BN$K.ID))], BN$S.cc[which(grepl(pattern = "PDGF", x = BN$S.cc))]))

idx1 <- which(BN$S.ID%in%allPDGF)
if(length(idx1) > 0){BN$S.ID[idx1] <- "PDGFR"}

idx2 <- which(BN$K.ID%in%allPDGF)
if(length(idx2)>0){BN$K.ID[idx2] <- "PDGFR"}

idx3 <- which(BN$S.cc%in%allPDGF)
if(length(idx3)>0){for(i in 1:length(idx3)){BN$S.cc[idx3[i]]<-paste0("PDGF_", strsplit(x = BN$S.cc[idx3[i]], split = "_")[[1]][2])}}

uniprot <- BN$S.AC[which(BN$S.ID=="PDGFR")[1]]

BN$S.AC[idx1] <- "PDGFR"

BN$K.AC[idx2] <- "PDGFR"

###
idx <- c()
for(i in 1:nrow(BN)){
  
  BN$S.AC[i] <- BN$S.ID[i]
  BN$K.AC[i] <- BN$K.ID[i]
  
  if(BN$S.ID[i]==BN$K.ID[i]){
    
    idx <- c(idx, i)
    
  }
  
}

BN <- BN[-idx, ]

# idx <- intersect(x = which(BN$K.ID=="PIK3C"), y = which(BN$S.ID=="AKT"))
# BN <- BN[-idx, ]

idx <- intersect(x = which(BN$K.ID=="PIK3C"), y = which(BN$S.ID=="MTOR")) # replace PIK3C->MTOR with PIK3C-PDPK1
BN <- BN[-idx, ]
# BN$S.AC[idx] <- "PDPK1"
# BN$S.ID[idx] <- "PDPK1"
# BN$S.cc[idx] <- "PDPK1_R1"

#
save(BN, file = "allD_MOUSE.RData")
