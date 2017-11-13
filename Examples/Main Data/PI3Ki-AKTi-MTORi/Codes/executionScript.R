# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("buildDataMatrix.R")
source("ilpFunctions.R")
source("ilpFunctions2.R")

# load(file='omnipath_allD.Rdata')
load(file='allD_noCSK_filt.RData')
load(file='dataObjects_PHONEMeS.RData')

#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM.res.noFC, IDmap=GMM.res.ID, resFC=GMM.res)

conditions <- list(c("AKT1 - Control", "AKT2 - Control"), c("CAMK1 - Control", "CAMK2 - Control"),
                    c("EGFR1 - Control", "EGFR2 - Control"), c("ERK1 - Control", "ERK2 - Control"),
                    c("MEK1 - Control", "MEK2 - Control"), c("MTOR1 - Control", "MTOR2 - Control"),
                    c("P70S6K1 - Control", "P70S6K2 - Control"), c("PI3K2 - Control", "PI3K2 - Control"),
                    c("PKC1 - Control", "PKC2 - Control"), c("ROCK1 - Control", "ROCK2 - Control"))

names(conditions) <- c("AKT1_HUMAN", "KCC2D_HUMAN", "EGFR_HUMAN", "MK01_HUMAN",
                        "MP2K1_HUMAN", "MTOR_HUMAN",  "KS6B1_HUMAN", "PK3CA_HUMAN",
                        "KPCA_HUMAN", "ROCK1_HUMAN")

#Choose the drug targets
resultsSIF <- matrix(, nrow = 1, ncol = 3)
targets.P<-list(cond1=c("MTOR_HUMAN", "AKT1_HUMAN", "PK3CA_HUMAN"))
for(ii in 1:length(targets.P[[1]])){
  
  targets <- list(cond1=c(targets.P[[1]][ii]))
  
  #Choose the drug treatments matching to the drug targets
  #and match to what is present in the background network
  # data.P<-dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=list(cond1=c("MTOR1 - Control", "MTOR2 - Control", "AKT1 - Control", "AKT2 - Control")))
  # experiments <- list(cond1=c("MTOR1 - Control", "MTOR2 - Control", "AKT1 - Control", "AKT2 - Control"))
  # data.P<-dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=list(cond1=rownames(GMM.res[[1]])[c(2,3,12,13,16,17)]))
  # experiments <- list(cond1=rownames(GMM.res[[1]])[c(2,3,12,13,16,17)])
  data.P <- dataBycond(dataGMM, bg, scaled = TRUE, rowBycond = conditions[which(names(conditions)==targets[[1]])])
  experiments <- conditions[which(names(conditions)==targets[[1]])]
  
  show(data.P)
  
  speciesP(data.P)
  #Create the PKN list that will be used for optimisation
  pknList<-buildNw(data.On=data.P, targets.On=targets, bg=bg,nK="no")
  idx <- which(pknList@interactions$S.ID==pknList@interactions$K.ID)
  rem <- c()
  if(length(idx) > 0){
    
    for(i in 1:length(idx)){
      
      # kinase <- pknList@interactions[idx[i], ]$K.ID
      # substrate <- pknList@interactions[idx[i], ]$S.cc
      
      # rem <- c(rem, intersect(which(pknList@interactions$K.ID==pknList@interactions[idx[1], ]$S.cc), which(pknList@interactions$S.cc==pknList@interactions[idx[1], ]$K.ID)))
      
      pknList@interactions <- pknList@interactions[-intersect(which(pknList@interactions$K.ID==pknList@interactions[idx[i], ]$S.cc), which(pknList@interactions$S.cc==pknList@interactions[idx[i], ]$K.ID)), ]
      
    }
    
  }
  
  
  show(pknList)
  
  # dataGMM@IDmap <- dataGMM@IDmap[-which(duplicated(dataGMM@IDmap$S.cc)), ]
  # for(i in 1:length(which(duplicated(dataGMM@IDmap$S.cc)))){
  #   
  #   dataGMM@res[-which(names(dataGMM@res)==dataGMM@IDmap$dataID[which(duplicated(dataGMM@IDmap$S.cc))[i]])]
  #   dataGMM@resFC[-which(names(dataGMM@res)==dataGMM@IDmap$dataID[which(duplicated(dataGMM@IDmap$S.cc))[i]])]
  #   
  # }
  
  
  ################################################################################
  #                         First step optimization                              #
  ################################################################################
  # Build the matrix wth the necessary data for all the species in the prior knowledge
  dataMatrix <- buildDataMatrix(dataGMM = dataGMM, pknList = pknList, targets = targets, experiments = experiments)
  
  # SIF file for the prior knowledge interactions
  sif <- createSIF(pknList = pknList)
  
  # Creating lists containing all our variables
  binary_x <- create_binary_variables_for_x_vector(dataMatrix = dataMatrix)
  binary_y <- create_binary_variables_for_y_vector(pknList = pknList)
  binaries <- create_binaries(binaries_x = binary_x, binaries_y = binary_y)
  
  # Writing the objective function
  oF <- write_objective_function(dataMatrix = dataMatrix, binaries = binaries)
  
  # Writing the bounds and also all the vvariables
  bounds <- write_boundaries(binaries = binaries)
  
  # Writing equality constraints
  eC <- write_equality_constraints(dataMatrix = dataMatrix, binaries = binaries)
  
  # Writing Constraint - 1
  c1 <- write_constraints_1(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Writing Constraint - 2
  c2 <- write_constraints_2(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Writing Constraint - 3
  c3 <- write_constraints_3(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Writing Constraint - 4
  c4 <- write_constraints_4(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Writeing Constraint - 5
  c5 <- write_constraints_5(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Putting all constraints together in one file
  allC <- all_constraints(equalityConstraints = eC, constraints1 = c1, constraints2 = c2, constraints3 = c3, constraints4 = c4, constraints5 = c5)
  
  # write(bounds, file = "bounds.txt")
  # write(binaries[[1]], file = "Integers.txt")
  # write(allC, file = "allConstraints.txt")
  
  # write the .lp file
  data = "testFile.lp"
  write("enter Problem", data)
  write("", data, append = TRUE)
  write("Minimize", data, append = TRUE)
  write(oF, data, append = TRUE)
  write("Subject To", data, append = TRUE)
  write(allC, data, append = TRUE)
  write("Bounds", data, append = TRUE)
  write(bounds, data, append = TRUE)
  write("Binaries", data, append = TRUE)
  write(binaries[[1]], data, append = TRUE)
  write("End", data, append = TRUE)
  
  system(paste0(getwd(), "/cplex -f cplexCommand.txt"))
  
  
  # Read the results from the CPLEX and do the necessary processing of the model
  library(XML)
  resultsSIF1 <- readOutSIF(cplexSolutionFileName = "results1.txt", binaries = binaries)
  write.table(resultsSIF1, file = "resultsSIF1.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ################################################################################
  #                         Second step optimization                             #
  ################################################################################
  dataMatrix2 <- buildDataMatrix2(dataMatrix = dataMatrix, sif = resultsSIF1, targets = targets, experiments = experiments)
  binX <- create_binary_variables_for_x_vector_step_2(dataMatrix = dataMatrix2)
  binY <- create_binary_variables_for_y_vector_step_2(sif = resultsSIF1)
  binaries2 <- create_binaries(binaries_x = binX, binaries_y = binY)
  oF2 <- write_objective_function_step_2Edges(binaries = binaries2, dataMatrix = dataMatrix2)
  bounds2 <- write_boundaries_step_2(binaries = binaries2)
  c1 <- write_fixed_nodes_constraints(dataMatrix = dataMatrix2, binaries = binaries2)
  c2 <- write_one_edge_constraint(sif = resultsSIF1, dataMatrix = dataMatrix2, binaries = binaries2)
  c3 <- write_intermediate_node_constraints_in(sif = resultsSIF1, dataMatrix = dataMatrix2, binaries = binaries2)
  c4 <- write_intermediate_node_constraints_out(sif = resultsSIF1, dataMatrix = dataMatrix2, binaries = binaries2)
  c5 <- write_reaction_constraints(dataMatrix = dataMatrix2, sif = resultsSIF1, binaries = binaries2)
  c6 <- write_loop_constraints(resultsSIF1 = resultsSIF1, binaries = binaries2)
  allC <- all_constraints_2(constraints1 = c1, constraints2 = c2, constraints3 = c3, constraints4 = c4, constraints5 = c5, constraints6 = c6)
  
  # write the .lp file
  data = "testFile2.lp"
  write("enter Problem", data)
  write("", data, append = TRUE)
  write("Minimize", data, append = TRUE)
  write(oF2, data, append = TRUE)
  write("Subject To", data, append = TRUE)
  write(allC, data, append = TRUE)
  write("Bounds", data, append = TRUE)
  write(bounds2, data, append = TRUE)
  write("Binaries", data, append = TRUE)
  write(binaries2[[1]], data, append = TRUE)
  write("End", data, append = TRUE)
  
  system(paste0(getwd(), "/cplex -f cplexCommand2.txt"))
  
  library(CellNOptR)
  resultsSIF2 <- readOutSIF(cplexSolutionFileName = "results2.txt", binaries = binaries2)
  resultsSIF <- rbind(resultsSIF, resultsSIF2)

  file.remove("cplex.log")
  # file.remove("clone1.log")
  file.remove("results1.txt")
  file.remove("results2.txt")
  file.remove("testFile.lp")
  file.remove("testFile2.lp")
  
}

write.table(unique(resultsSIF[-1, ]), file = "resultsSIF.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
