# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)
library(dplyr)

# Call the PHONEMeS functions
source("../../../Public/buildDataMatrix.R")
source("../../../Public/ilpFunctions.R")
source("../../../Public/buildDataObject.R")
source("../../../Public/build_Nw.R")
source("../../../Public/build_PKN.R")

# Loading BN and Data
load(file = "allD_noCSK_filt.RData")
load(file = "data_Hsu_PHONEMeS.RData")

# Building the PHONEMeS inpputs
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)

#Choose the drug targets
targets.P <- list(cond1=c("MTOR_HUMAN"), cond2=c("MTOR_HUMAN"))

# Choose the drug treatments matching to the drug targets and match to what is present in the background network
conditions <- list(cond1=c("rapa - ins"), cond2=c("tor - ins"))

# Select experimental condition
experiments <- c(2) # for tor -- experiment <- c(1) for rapa

# Run PHONEMeS with multiple solutions
source("../../../Public/runPHONEMeS.R")
resultsMulti <- runPHONEMeS(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg, solver = "cplex", nSolutions = 100, nK = "no")
write.table(x = resultsMulti, file = "tor_sif_cplex.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

source("../../../Public/runPHONEMeS_Downsampling.R")
resultsSampling <- runPHONEMeS_Downsampling(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg, nIter = 100, solver = "cplex", nK = "no")
write.table(x = resultsSampling, file = "tor_sif_downsampling.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = resultsSampling[which(resultsSampling[, 2]>=10), ], file = "tor_sif_downsampling_cutoff10.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Assign nodes attributes for nicer visualization in cytoscape
source("../../../Public/assignAttributes.R")
nodesAttributes <- assignAttributes(sif = resultsSampling, dataGMM = dataGMM, targets = targets.P[experiments], writeAttr = TRUE)
