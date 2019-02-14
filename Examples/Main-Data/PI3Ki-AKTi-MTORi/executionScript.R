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

# Loading database and dat-object
load(file='allD_noCSK_filt.RData')
load(file='dataObjects_PHONEMeS.RData')

#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM.res.noFC, IDmap=GMM.res.ID, resFC=GMM.res)

conditions <- list(c("AKT1 - Control", "AKT2 - Control"), c("CAMK1 - Control", "CAMK2 - Control"),
                   c("EGFR1 - Control", "EGFR2 - Control"), c("ERK1 - Control", "ERK2 - Control"),
                   c("MEK1 - Control", "MEK2 - Control"), c("MTOR1 - Control", "MTOR2 - Control"),
                   c("P70S6K1 - Control", "P70S6K2 - Control"), c("PI3K1 - Control", "PI3K2 - Control"),
                   c("PKC1 - Control", "PKC2 - Control"), c("ROCK1 - Control", "ROCK2 - Control"))

names(conditions) <- c("AKT1_HUMAN", "KCC2D_HUMAN", "EGFR_HUMAN", "MK01_HUMAN",
                       "MP2K1_HUMAN", "MTOR_HUMAN",  "KS6B1_HUMAN", "PK3CA_HUMAN",
                       "KPCA_HUMAN", "ROCK1_HUMAN")

targets.P<-list(cond1=c("AKT1_HUMAN", "AKT2_HUMAN"), cond2=c("KCC2A_HUMAN", "KCC2B_HUMAN", "KCC2C_HUMAN", "KCC2D_HUMAN"), cond3=c("EGFR_HUMAN", "ERBB2_HUMAN"), 
                cond4=c("MK01_HUMAN", "MK03_HUMAN", "MK14_HUMAN"), cond5=c("MP2K1_HUMAN", "MP2K2_HUMAN"), cond6=c("MTOR_HUMAN"), cond7=c("KS6B1_HUMAN", "KS6B2_HUMAN"), 
                cond8=c("PK3CA_HUMAN", "PK3CD_HUMAN", "MTOR_HUMAN"), cond9=c("KPCA_HUMAN", "KPCB_HUMAN", "KPCG_HUMAN", "KPCE_HUMAN"), cond10=c("ROCK1_HUMAN", "ROCK2_HUMAN"))

# Select experimental condition
experiments <- c(1, 6, 8) # for PI3Ki-AKTi-MTORi case

# Running PHONEMeS - cplex
# Run PHONEMeS with multiple solutions from CPLEX
source("../../../Public/runPHONEMeS.R")
resultsMulti <- runPHONEMeS(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg, solver = "cplex", nSolutions = 100, nK = "no")
write.table(x = resultsMulti, file = "Cascade_sif_cplex.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Run PHONEMeS with multiple solutions from Downsampling
source("../../../Public/runPHONEMeS_Downsampling.R")
resultsSampling <- runPHONEMeS_Downsampling(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg, solver = "cplex", nIter = 100, nK = "no")
write.table(x = resultsSampling, file = "Cascade_sif_sampling.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = resultsSampling[which(resultsSampling[, 2]>10), ], file = "Cascade_sif_cascade_cutoff10.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Assign nodes attributes for nicer visualization in cytoscape
source("../../../Public/assignAttributes.R")
nodesAttributes <- assignAttributes(sif = resultsSampling, dataGMM = dataGMM, targets = targets.P[experiments], writeAttr = TRUE)
