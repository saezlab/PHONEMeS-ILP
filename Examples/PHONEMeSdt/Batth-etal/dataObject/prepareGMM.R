# Setting the parameters for the input
source("../../../../Public/buildInputs.R")
load(file = "ttop_list.RData")

# Highly regulated sites are considered as perturbed
idxFC = 2
idxPval = 6
idxID = 1
pThresh = rep(0.025, 6)
fcThresh = rep(2.5, 6)
namesConditions = c("3min_PDGF", "3min_IGF1", "3min_FGF2", "15min_PDGF", "15min_IGF1", "15min_FGF2")

inputs <- buildInputs(tableTopList = ttop_list, fcThresh = fcThresh, pThresh = pThresh, idxID = idxID, idxFC = idxFC, idxPval = idxPval, namesConditions = namesConditions)
GMM.ID <- inputs$IDmap
GMM <- inputs$res
GMM.wFC <- inputs$resFC

save(GMM.ID, file = "GMM.ID.RData")
save(GMM, file = "GMM.RData")
save(GMM.wFC, file = "GMM.wFC.RData")
