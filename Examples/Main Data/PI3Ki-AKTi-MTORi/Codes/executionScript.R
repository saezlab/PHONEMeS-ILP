#
#  This file is part of the CNO software
#
#  Copyright (c) 2018 - RWTH Aachen - JRC COMBINE
#
#  File author(s): E. Gjerga (enio.gjerga@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: https://saezlab.github.io/PHONEMeS/
#
##############################################################################
# $Id$

# Preparing & Executing PHONEMeS

# Load packages
library("BioNet")
library("igraph")
library("PHONEMeS")
library("hash")
library("dplyr")
library("tidyr")

# Load functions
source("buildDataMatrix.R")
source("ilpFunctions.R")
source("build_Nw.R")
source("build_PKN.R")
source("runPHONEMeS.R")

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

targets.P<-list(cond1=c("AKT1_HUMAN", "AKT2_HUMAN"), cond2=c(), cond3=c(), cond4=c(),
                cond5=c(), cond6=c("MTOR_HUMAN"), cond7=c(), cond8=c("PK3CA_HUMAN", "PK3CD_HUMAN", "MTOR_HUMAN"),
                cond9=c(), cond10=c())

# Running PHONEMeS - cplex
resultsSIF_cplex <- runPHONEMeS(targets.P=targets.P, 
                                conditions=conditions, 
                                dataGMM=dataGMM, 
                                experiments=c(1, 6, 8), 
                                bg=bg, 
                                nK="no", 
                                solver="cplex")

# Writing complete set of results
write.table(x = resultsSIF_cplex, file = "resultSIF_cplex.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Reducing result network to include the most likely interactions
source("reduceSIF.R")
reducedSIF <- reduceSIF(sif = resultsSIF_cplex, targets = targets.P, dataGMM = dataGMM, cutoff = 0.05)
write.table(x = reducedSIF, file = "reducedSIF_005.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Assigning nodes attributes for better visualization
source("assignAttributes.R")
nodesAttributes <- assignAttributes(sif = resultsSIF_cplex, dataGMM = dataGMM, targets = targets.P, writeAttr = TRUE)

# Running PHONEMeS - cbc
resultsSIF_cbc <- runPHONEMeS(targets.P=targets.P, 
                              conditions=conditions, 
                              dataGMM=dataGMM, 
                              experiments=c(6), 
                              bg=bg, 
                              nK="all", 
                              solver="cbc")

write.table(resultsSIF_cplex, file = "resultsSIF.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# compare cplex and cbc results
resultsSIF_cplex$solver <- "cplex"
resultsSIF_cbc$solver <- "cbc"
resultsSIF <- bind_rows(resultsSIF_cplex, resultsSIF_cbc)
tmp <- resultsSIF %>% group_by(Source, Interaction, Target) %>% summarise(n=n()) %>% ungroup()
