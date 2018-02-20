# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)
library(readr)

##
load("tableTopList.RData")
mappingTable <- read_csv("mappingTable.csv")
mappingTable <- mappingTable[, -1] # first column is redundant
source("buildInputs.R")

dataInput <- buildInputs(tableTopList = ttList, organism = "HUMAN", pThresh = 0.1, mappingTable = mappingTable)

##
load("allD_HUMAN_Omnipath.RData")
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=dataInput$res, IDmap=dataInput$IDmap, resFC=dataInput$resFC)

##
targets.P <- list(cond1=c("MTOR_HUMAN"))

conditions <- list(c("cond1"))

names(conditions) <- "MTOR - Control"

data.P <- dataBycond(dataGMM, bg, scaled = TRUE, rowBycond = conditions)
experiments <- conditions

show(data.P)

speciesP(data.P)

##
source("buildDataMatrix.R")
source("ilpFunctions.R")
source("PKN_list.R")
source("build_Nw_Inv.R")
source("build_PKN_Inv.R")
source("runPHONEMeS.R")
source("runPHONEMeS-UD.R")

res <- runPHONEMeS(dataInput = dataInput, bg = bg, targets.P = targets.P, data.P = data.P, experiments = experiments)
resUD <- runPHONEMeSUD(dataInput = dataInput, bg = bg, targets.P = targets.P, data.P = data.P, experiments = experiments)
