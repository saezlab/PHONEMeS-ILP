# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)
library(dplyr)

# Call the PHONEMeS functions
source("../../Public/buildDataMatrix.R")
source("../../Public/ilpFunctions.R")
source("../../Public/buildDataObject.R")
source("../../Public/build_Nw.R")
source("../../Public/build_PKN.R")
source("../../Public/build_Nw_Inv.R")
source("../../Public/build_PKN_Inv.R")
source("../../Public/PKN_list.R")
source("../../Public/assignAttributes.R")
source("../../Public/runPHONEMeS_UD.R")

load(paste0(getwd(), "/Background-Network/allD.RData"))
load(paste0(getwd(), "/Data/GMM.ID.RData"))
load(paste0(getwd(), "/Data/GMM.RData"))
load(paste0(getwd(), "/Data/GMM.wFC.RData"))
load(paste0(getwd(), "/Data/kinaseList.RData"))

#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)

conditions <- list(c("cancer_vs_healthy"), c("metastasis_vs_healthy"), c("metastasis_vs_cancer"), c("metastasis_resolution"))
names(conditions) <- c("cond1", "cond2", "cond3", "cond4")

targets.P<-list(cond1=kinaseList[[1]], cond2=kinaseList[[2]], cond3=kinaseList[[3]], cond4=kinaseList[[4]])

# Select experimental condition
experiments <- c(1) # cancer_vs_healthy
sif <- runPHONEMeS_UD(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg)
write.table(x = sif[[1]], file = "down_cancer_vs_healthy.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
if(!is.null(sif[[2]])){write.table(x = sif[[2]], file = "up_cancer_vs_healthy.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)}
write.table(x = sif[[3]], file = "comb_cancer_vs_healthy.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
nodesAttributes <- assignAttributes(sif = sif[[3]], dataGMM = dataGMM, targets = targets.P[[experiments]], writeAttr = TRUE, fileName = "attributes_cancer_vs_healthy.txt")
rm(sif)

experiments <- c(2) # metastasis_vs_healthy
sif <- runPHONEMeS_UD(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg)
write.table(x = sif[[1]], file = "down_metastasis_vs_healthy.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
if(!is.null(sif[[2]])){write.table(x = sif[[2]], file = "up_metastasis_vs_healthy.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)}
write.table(x = sif[[3]], file = "comb_metastasis_vs_healthy.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
nodesAttributes <- assignAttributes(sif = sif[[3]], dataGMM = dataGMM, targets = targets.P[[experiments]], writeAttr = TRUE, fileName = "attributes_metastasis_vs_healthy.txt")
rm(sif)

experiments <- c(3) # metastasis_vs_cancer
sif <- runPHONEMeS_UD(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg)
write.table(x = sif[[1]], file = "down_metastasis_vs_cancer.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
if(!is.null(sif[[2]])){write.table(x = sif[[2]], file = "up_metastasis_vs_cancer.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)}
write.table(x = sif[[3]], file = "comb_metastasis_vs_cancer.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
nodesAttributes <- assignAttributes(sif = sif[[3]], dataGMM = dataGMM, targets = targets.P[[experiments]], writeAttr = TRUE, fileName = "attributes_metastasis_vs_cancer.txt")
rm(sif)

experiments <- c(4)
sif <- runPHONEMeS_UD(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg)
write.table(x = sif[[1]], file = "down_metastasis_resolution.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
if(!is.null(sif[[2]])){write.table(x = sif[[2]], file = "up_metastasis_resolution.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)}
write.table(x = sif[[3]], file = "comb_metastasis_resolution.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
nodesAttributes <- assignAttributes(sif = sif[[3]], dataGMM = dataGMM, targets = targets.P[[experiments]], writeAttr = TRUE, fileName = "attributes_metastasis_resolution.txt")
rm(sif)
