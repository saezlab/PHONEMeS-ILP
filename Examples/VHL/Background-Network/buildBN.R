# Building the background network from online resources of Omnipath (k/p-s only for this case)

library(readr)

url <- paste0(
  'http://omnipathdb.org/ptms?',
  'fields=sources,references&genesymbols=1'
)

download_omnipath <- function(){
  
  read.table(url, sep = '\t', header = TRUE)
  
}

omnipath_ptm <- download_omnipath()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]

##
BN <- matrix(data = , nrow = nrow(omnipath_ptm), ncol = 8)
colnames(BN) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

BN[, 1] <- as.character(omnipath_ptm$substrate_genesymbol)
BN[, 2] <- as.character(omnipath_ptm$substrate_genesymbol)
BN[, 3] <- as.character(omnipath_ptm$enzyme_genesymbol)
BN[, 4] <- as.character(omnipath_ptm$enzyme_genesymbol)
BN[, 5] <- as.character(omnipath_ptm$residue_type)
BN[, 6] <- as.character(omnipath_ptm$residue_offset)
BN[, 7] <- paste0("e", 1:nrow(omnipath_ptm))
BN[, 8] <- paste0(as.character(omnipath_ptm$substrate_genesymbol), "_", as.character(omnipath_ptm$residue_type), as.character(omnipath_ptm$residue_offset))

allD <- as.data.frame(BN)
allD$S.AC <- as.character(allD$S.AC)
allD$S.ID <- as.character(allD$S.ID)
allD$K.AC <- as.character(allD$K.AC)
allD$K.ID <- as.character(allD$K.ID)
allD$res <- as.character(allD$res)
allD$pos <- as.character(allD$pos)
allD$SID <- as.character(allD$SID)
allD$S.cc <- as.character(allD$S.cc)

save(allD, file = "allD.RData")
