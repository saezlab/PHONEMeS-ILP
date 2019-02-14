library(readr)

KSN_mouse_uniprot <- read_delim("KSN_mouse_uniprot", "\t", escape_double = FALSE, trim_ws = TRUE)

omnipath_mouse <- read_delim("omnipath_mouse", "\t", escape_double = FALSE, trim_ws = TRUE)

mappingTable <- matrix(data = , nrow = 2*nrow(omnipath_mouse), ncol = 2)
colnames(mappingTable) <- c("Uniprot", "Gene")
mappingTable[, 1] <- c(omnipath_mouse$UniProt_A, omnipath_mouse$UniProt_B)
mappingTable[, 2] <- c(omnipath_mouse$GeneSymbol_A, omnipath_mouse$GeneSymbol_B)
mappingTable <- unique(mappingTable)

speciesKSN <- unique(c(KSN_mouse_uniprot$enzyme, KSN_mouse_uniprot$substrate))
idx <- which(speciesKSN%in%mappingTable[, 1])
speciesKSN <- speciesKSN[-idx]

write.table(x = as.matrix(speciesKSN), file = "toUniprot.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE) # map uniprot identifiers to gene ID's and save the mapping table locally

toBind <- matrix(data = , nrow = length(speciesKSN), ncol = 2)
colnames(toBind) <- c("Uniprot", "Gene")
toBind[, 1] <- speciesKSN
uniprot <- read_delim("uniprot", "\t", escape_double = FALSE, trim_ws = TRUE)
for(i in 1:nrow(toBind)){
  
  if(length(which(uniprot$Entry==toBind[i, 1])) > 0){
    
    toBind[i, 2] <- strsplit(x = uniprot$`Gene names`[which(uniprot$Entry==toBind[i, 1])], split = " ", fixed = TRUE)[[1]][1]
    
  } else {
    
    toBind[i, 2] <- toBind[i, 1]
  }

}

mappingTable <- rbind(mappingTable, toBind)

mappingTable[which(is.na(mappingTable[, 2])), 2] <- mappingTable[which(is.na(mappingTable[, 2])), 1]

##
# Adding kinase/phosphatase-substrate interactions
BN <- matrix(data = , nrow = nrow(KSN_mouse_uniprot), ncol = 8)
colnames(BN) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

cnt <- 1
for(i in 1:nrow(KSN_mouse_uniprot)){
  
  BN[i, 1] <- KSN_mouse_uniprot$substrate[i]
  if(length(which(mappingTable[, 1]==KSN_mouse_uniprot$substrate[i])) > 0){
    BN[i, 2] <- toupper(mappingTable[which(mappingTable[, 1]==KSN_mouse_uniprot$substrate[i]), 2])
  } else { BN[i, 2] <- BN[i, 1] }
  BN[i, 3] <- KSN_mouse_uniprot$enzyme[i]
  if(length(which(mappingTable[, 1]==KSN_mouse_uniprot$enzyme[i])) > 0){
    BN[i, 4] <- toupper(mappingTable[which(mappingTable[, 1]==KSN_mouse_uniprot$enzyme[i]), 2])
  } else { BN[i, 4] <- BN[i, 3] }
  BN[i, 5] <- KSN_mouse_uniprot$residue[i]
  BN[i, 6] <- KSN_mouse_uniprot$offset[i]
  BN[i, 7] <- paste0("e", cnt)
  cnt <- cnt + 1
  BN[i, 8] <- paste0(BN[i, 2], "_", BN[i, 5], BN[i, 6])
  
}

BN1 <- BN

# ##
# # Adding the signed and directed interactions from uniprot
# BN <- matrix(data = , nrow = nrow(omnipath_mouse), ncol = 8)
# colnames(BN) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")
# 
# for(i in 1:nrow(omnipath_mouse)){
#   
#   BN[i, 1] <- omnipath_mouse$UniProt_B[i]
#   BN[i, 2] <- toupper(omnipath_mouse$GeneSymbol_B[i])
#   BN[i, 3] <- omnipath_mouse$UniProt_A[i]
#   BN[i, 4] <- toupper(omnipath_mouse$GeneSymbol_A[i])
#   BN[i, 5] <- "R"
#   BN[i, 6] <- "1"
#   BN[i, 7] <- paste0("e", cnt)
#   cnt <- cnt + 1
#   BN[i, 8] <- paste0(BN[i, 2], "_R1")
#   
# }
# 
# BN2 <- BN
# 
# ##
# BN <- rbind(BN1, BN2)

##
# add igfr1 from string interactions
string_interactions_igf1r <- read_delim("string_interactions_igf1r.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
interactions <- string_interactions_igf1r[which(string_interactions_igf1r$node2=="Igf1r"), ]
interactions$`#node1`[which(interactions$`#node1`=="Hras1")] <- "Hras"
for(i in 1:nrow(interactions)){
  
  if((i%in%c(2, 4)) == FALSE){
    
    cc1 <- as.character(mappingTable[which(mappingTable[, 2]==interactions$`#node1`[i]), 1])
    cc2 <- toupper(interactions$`#node1`[i])
    cc3 <- as.character(mappingTable[which(mappingTable[, 2]==interactions$node2[i]), 1])
    cc4 <- toupper(interactions$node2[i])
    cc5 <- "R"
    cc6 <- "1"
    cc7 <- paste0("e", cnt)
    cnt <- cnt + 1
    cc8 <- paste0(cc2, "_R1")
    
    toBind <- as.data.frame(t(as.matrix(c(cc1, cc2, cc3, cc4, cc5, cc6, cc7, cc8))))
    colnames(toBind) <- colnames(BN)
    
    BN <- rbind(BN, toBind)
    
  }
  
}

interactions <- string_interactions_igf1r[which(string_interactions_igf1r$`#node1`=="Igf1r"), ]
for(i in 1:nrow(interactions)){
  
  cc1 <- as.character(mappingTable[which(mappingTable[, 2]==interactions$node2[i]), 1])
  cc2 <- toupper(interactions$node2[i])
  cc3 <- as.character(mappingTable[which(mappingTable[, 2]==interactions$`#node1`[i]), 1])
  cc4 <- toupper(interactions$`#node1`[i])
  cc5 <- "R"
  cc6 <- "1"
  cc7 <- paste0("e", cnt)
  cnt <- cnt + 1
  cc8 <- paste0(cc2, "_R1")
  
  toBind <- as.data.frame(t(as.matrix(c(cc1, cc2, cc3, cc4, cc5, cc6, cc7, cc8))))
  colnames(toBind) <- colnames(BN)
  
  BN <- rbind(BN, toBind)
  
}

##
# add fgfr1 from string interactions
string_interactions_fgfr1 <- read_delim("string_interactions_fgfr1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
interactions <- string_interactions_fgfr1[which(string_interactions_fgfr1$`#node1`=="Fgfr1"), ]
interactions$node2[which(interactions$node2=="Hras1")] <- "Hras"
for(i in 1:nrow(interactions)){
  
  if(i!=1){
    
    cc1 <- as.character(mappingTable[which(mappingTable[, 2]==interactions$node2[i]), 1])
    cc2 <- toupper(interactions$node2[i])
    cc3 <- as.character(mappingTable[which(mappingTable[, 2]==interactions$`#node1`[i]), 1])
    cc4 <- toupper(interactions$`#node1`[i])
    cc5 <- "R"
    cc6 <- "1"
    cc7 <- paste0("e", cnt)
    cnt <- cnt + 1
    cc8 <- paste0(cc2, "_R1")
    
    toBind <- as.data.frame(t(as.matrix(c(cc1, cc2, cc3, cc4, cc5, cc6, cc7, cc8))))
    colnames(toBind) <- colnames(BN)
    
    BN <- rbind(BN, toBind)
    
  }
  
}

##
# add pdgfra from string interactions
string_interactions_pdgfra <- read_delim("string_interactions_pdgfra.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
interactions <- string_interactions_pdgfra[which(string_interactions_pdgfra$node2=="Pdgfra"), ]
interactions$`#node1`[which(interactions$`#node1`=="Hras1")] <- "Hras"
for(i in 1:nrow(interactions)){
  
  if((i%in%c(1, 2, 6, 10)) == FALSE){
    
    cc1 <- as.character(mappingTable[which(mappingTable[, 2]==interactions$`#node1`[i]), 1])
    cc2 <- toupper(interactions$`#node1`[i])
    cc3 <- as.character(mappingTable[which(mappingTable[, 2]==interactions$node2[i]), 1])
    cc4 <- toupper(interactions$node2[i])
    cc5 <- "R"
    cc6 <- "1"
    cc7 <- paste0("e", cnt)
    cnt <- cnt + 1
    cc8 <- paste0(cc2, "_R1")
    
    toBind <- as.data.frame(t(as.matrix(c(cc1, cc2, cc3, cc4, cc5, cc6, cc7, cc8))))
    colnames(toBind) <- colnames(BN)
    
    BN <- rbind(BN, toBind)
    
  }
  
}

##
# add pdgfrb from string interactions
string_interactions_pdgfrb <- read_delim("string_interactions_pdgfrb.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
interactions <- string_interactions_pdgfrb[which(string_interactions_pdgfrb$node2=="Pdgfrb"), ]
for(i in 1:nrow(interactions)){
  
  if((i%in%c(6, 7, 8)) == FALSE){
    
    cc1 <- as.character(mappingTable[which(mappingTable[, 2]==interactions$`#node1`[i]), 1])
    cc2 <- toupper(interactions$`#node1`[i])
    cc3 <- as.character(mappingTable[which(mappingTable[, 2]==interactions$node2[i]), 1])
    cc4 <- toupper(interactions$node2[i])
    cc5 <- "R"
    cc6 <- "1"
    cc7 <- paste0("e", cnt)
    cnt <- cnt + 1
    cc8 <- paste0(cc2, "_R1")
    
    toBind <- as.data.frame(t(as.matrix(c(cc1, cc2, cc3, cc4, cc5, cc6, cc7, cc8))))
    colnames(toBind) <- colnames(BN)
    
    BN <- rbind(BN, toBind)
    
  }
  
}

##
# Add PI3KCA -> AKT1 since it is not present as ppi
cc1 <- "P31750"
cc2 <- "AKT1"
cc3 <- "P42337"
cc4 <- "PIK3CA"
cc5 <- "R"
cc6 <- "1"
cc7 <- paste0("e", cnt)
cc8 <- "AKT1_R1"

toBind <- as.data.frame(t(as.matrix(c(cc1, cc2, cc3, cc4, cc5, cc6, cc7, cc8))))
colnames(toBind) <- colnames(BN)

BN <- rbind(BN, toBind)

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
BN$S.ID[which(BN$S.ID=="PDGFRB")] <- "PDGFR"
BN$S.ID[which(BN$S.ID=="PDGFRA")] <- "PDGFR"

BN$K.ID[which(BN$K.ID=="PDGFRB")] <- "PDGFR"
BN$K.ID[which(BN$K.ID=="PDGFRA")] <- "PDGFR"

BN$S.cc[which(grepl(pattern = "PDGFR", x = BN$S.cc))] <- gsub(pattern = "PDGFRA", replacement = "PDGFR", fixed = TRUE, x = BN$S.cc[which(grepl(pattern = "PDGFR", x = BN$S.cc))])
BN$S.cc[which(grepl(pattern = "PDGFR", x = BN$S.cc))] <- gsub(pattern = "PDGFRB", replacement = "PDGFR", fixed = TRUE, x = BN$S.cc[which(grepl(pattern = "PDGFR", x = BN$S.cc))])

BN <- BN[-c(12260, 12257), ]

##
idx <- intersect(x = which(BN$S.ID=="IGF1R"), y = which(BN$K.ID=="IGF1R"))
BN <- BN[-idx, ]

idx <- intersect(x = which(BN$S.ID=="FGFR1"), y = which(BN$K.ID=="FGFR1"))
BN <- BN[-idx, ]

idx <- intersect(x = which(BN$S.ID=="PDGFR"), y = which(BN$K.ID=="PDGFR"))
BN <- BN[-idx, ]

save(BN, file = "allD_MOUSE.RData")
