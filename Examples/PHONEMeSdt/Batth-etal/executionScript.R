#  Copyright (c) 2018 - RWTH Aachen University
#
#  File author(s): Enio Gjerga
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  email: enio.gjerga@gmai.com
#
##############################################################################
# 13:22 23/03/2018
# This script is used to produce the PHONEMeS Outputs by calling the functions
# which write the ILP formulation

#load the packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(CellNOptR)
library(tidyr)
library(dplyr)

# Call the PHONEMeS functions
source("../../../Public/buildDataMatrix.R")
source("../../../Public/ilpFunctions.R")
source("../../../Public/buildDataObject.R")
source("../../../Public/build_Nw.R")
source("../../../Public/build_PKN.R")
source("groupResidues.R")

load(paste0(getwd(), "/dataObject/GMM.ID.RData"))
load(paste0(getwd(), "/dataObject/GMM.RData"))
load(paste0(getwd(), "/dataObject/GMM.wFC.RData"))
load(paste0(getwd(), "/Background-Network/allD_MOUSE.RData"))

allD <- BN

#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)

conditions <- list(c("3min_PDGF"), c("3min_IGF1"), c("3min_FGF2"), c("15min_PDGF"), c("15min_IGF1"), c("15min_FGF2"))

names(conditions) <- c("cond1_3min", "cond2_3min", "cond3_3min", "cnod1_15min", "cond2_15min", "cond3_15min")

targets.P <- list(cond1_3min=c("PDGFR"), cond2_3min=c("IGF1R"), cond3_3min=c("FGFR2"), cond1_15min=c("PDGFR"), cond2_15min=c("IGF1R"), cond3_15min=c("FGFR2"))

experiments <- list(tp1=c(1), tp2=c(4)) # pdgfr only
# experiments <- list(tp1=c(1, 2, 3), tp2=c(4, 5, 6)) # if you want to consider all the three conditions for each time-point then do: experiments <- list(tp1=c(1, 2, 3), tp2=c(4, 5, 6))

#Generating the networks
source("../../../Public/runPHONEMeS_dt.R")
tpSIF <- runPHONEMeS_dt(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg, nIter = 100)
write.table(x = tpSIF, file = "pdgfr_tp_analysis.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=5), ], file = "pdgfr_tp_analysis_cutoff_5.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=10), ], file = "pdgfr_tp_analysis_cutoff_10.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=20), ], file = "pdgfr_tp_analysis_cutoff_20.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=50), ], file = "pdgfr_tp_analysis_cutoff_50.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=60), ], file = "pdgfr_tp_analysis_cutoff_60.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=70), ], file = "pdgfr_tp_analysis_cutoff_70.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=80), ], file = "pdgfr_tp_analysis_cutoff_80.txt", quote = FALSE, sep = "\t", row.names = FALSE)

source("../../../Public/assignAttributes_tp.R")
nodesAttributes <- assignAttributes_tp(sif = tpSIF, dataGMM = dataGMM, targets = targets.P[[1]], writeAttr = TRUE)

