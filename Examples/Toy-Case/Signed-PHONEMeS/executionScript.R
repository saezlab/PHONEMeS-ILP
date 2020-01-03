# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)
library(dplyr)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../../../Public/buildDataMatrix.R")
source("../../../Public/ilpFunctions.R")
source("../../../Public/buildDataObject.R")
source("../../../Public/build_Nw.R")
source("../../../Public/build_PKN.R")

load(file = "allD_Toy.RData")
load(file = "dataObjects_Toy.RData")

allD <- allD_Toy
GMM.res.noFC <- GMM.Toy
GMM.res <- GMM.Toy.wFC
GMM.res.ID <- GMM.Toy.ID

#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM.res.noFC, IDmap=GMM.res.ID, resFC=GMM.res)

conditions <- list()
for(i in 1:nrow(GMM.res[[1]])){
  
  conditions[[length(conditions)+1]] <- rownames(GMM.res[[1]])[i]
  
}

names(conditions) <- rownames(GMM.res[[i]])

targets.P<-list(cond1=c("S"), cond2=c("I"), cond3=c("S", "I"))

targets.S <- list(cond1=c("S"), cond2=c(""), cond3=c("S"))
targets.I <- list(cond1=c(""), cond2=c("I"), cond3=c("I"))

# targets.P <- targets.P[[1]]

source("../../../Public/runPHONEMeS_Signed.R")
resultsSIF_case1 <- runPHONEMeS_Signed(targets.S=targets.S,
                                targets.I=targets.I,       
                                conditions=conditions, 
                                dataGMM=dataGMM, 
                                experiments=c(1, 2, 3), 
                                bg=bg, 
                                solver="cplex")

write.table(x = resultsSIF_case1, file = "case1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
