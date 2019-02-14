source("../../../Public/buildInputs.R")

load(file = "../Preprocessing/input_PHONEMES_UD.Rdata")

ttop_list <- list()
for(ii in 1:length(for_enio)){
  
  ttop_list[[length(ttop_list)+1]] <- for_enio[[ii]]$ttop
  
}

# Highly regulated sites are considered as perturbed
idxFC = 2
idxPval = 6
idxID = 1
pThresh = c(0.025, 0.025, 0.05, 0.05)
fcThresh = c(3, 3, 2, 2)
namesConditions <- c("cancer_vs_healthy", "metastasis_vs_healthy", "metastasis_vs_cancer", "metastasis_resolution")

dataInput <- buildInputs(tableTopList = ttop_list, fcThresh = fcThresh, pThresh = pThresh, idxID = idxID, idxFC = idxFC, idxPval = idxPval, namesConditions = namesConditions)

GMM.ID <- dataInput$IDmap
GMM <- dataInput$res
GMM.wFC <- dataInput$resFC

save(GMM.ID, file = "GMM.ID.RData")
save(GMM, file = "GMM.RData")
save(GMM.wFC, file = "GMM.wFC.RData")

#
kinaseList <- list()
for(ii in 1:length(for_enio)){
  kinaseList[[length(kinaseList)+1]] <- for_enio[[ii]]$kinases
}
save(kinaseList, file = "kinaseList.RData")
