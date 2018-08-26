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

start.time <- Sys.time()
#load the packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(CellNOptR)

load("dataGMM.RData")
load("allD_MOUSE.RData")

allD <- BN
for(i in 1:length(GMM)){
  
  GMM[[i]] <- GMM[[i]][c(1, 4), ]
  GMM.wFC[[i]] <- GMM.wFC[[i]][c(1, 4), ]
  
}

source("build_Nw.R")
source("build_PKN.R")
source("buildDataMatrix.R")
source("ilpFunctions.R")
source("buildNetwork.R")

bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)

#Choose the drug targets
targets.P<-list(cond1=c("PDGFR"))
targets <- targets.P

resList <- list()

iterSum <- 0

GMM.all <- GMM
GMM.wFC.all <- GMM.wFC

SS <- list()

while(iterSum < 100){
  
  ss <- unique(sample(x = 1:length(GMM.all), replace = TRUE))
  SS[[length(SS)+1]] <- ss
  GMM <- list()
  GMM.wFC <- list()
  
  for(kk in 1:length(ss)){
    
    GMM[[length(GMM)+1]] <- GMM.all[[ss[kk]]]
    GMM.wFC[[length(GMM.wFC)+1]] <- GMM.wFC.all[[ss[kk]]]
    
  }
  
  names(GMM) <- names(GMM.all)[ss]
  names(GMM.wFC) <- names(GMM.all)[ss]
  
  bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
  dataGMM<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)
  
  resListSep <- list()
  
  for(i in 1:nrow(GMM[[1]])){
    
    if(i==1){
      
      delayedAssign("do.next", {next})
      
      experiments <- list(cond1=c(rownames(GMM[[1]])[i]))
      res <- tryCatch({dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=experiments)},
                      warning = function(w) {
                        force(do.next)
                      }, error = function(e) {
                        force(do.next)
                      })
      
      data.P<-dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=experiments)
      show(data.P)
      
      speciesP(data.P)
      
      res = tryCatch({build_Nw(data.On=data.P, targets.On=targets.P, bg=bg, nK = "yes")},
                     warning = function(w) {
                       force(do.next)
                     }, error = function(e) {
                       force(do.next)
                     })
      
      pknList <- build_Nw(data.On=data.P, targets.On=targets.P, bg=bg, nK = "yes")
      # cascade=list(c1=c("PDPK1_HUMAN"), c2=c("AKT_HUMAN"))
      # pknList <- addCascade(cascade = cascade, allD = allD, pknList = pknList)
      
      show(pknList)
      
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
      bounds <- write_boundaries(binaries = binaries, pknList = pknList, M = 100)
      
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
      
      # Writeing Constraint - 6
      c6 <- write_self_activating_constraints(pknList = pknList, binaries = binaries, dataMatrix = dataMatrix, M = 100)
      # c7 <- cascadeConstraints(pknList = pknList, cascade = cascade, binaries = binaries)
      # c6 <- c(c6, c7)
      
      # Putting all constraints together in one file
      allC <- all_constraints(equalityConstraints = eC, constraints1 = c1, constraints2 = c2, constraints3 = c3, constraints4 = c4, constraints5 = c5, constraints6 = c6)
      
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
      # write.table(resultsSIF1, file = "resultsSIF1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
      temp <- unique(resultsSIF1)
      resultsSIF1 <- removeRedundantNodes(resultsSIF1 = resultsSIF1)
      resultsSIF <- unique(resultsSIF1)
      colnames(resultsSIF) <- c("Source", "Interaction", "Target")
      write.table(resultsSIF, file = paste0("resultSIF_tp", i, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
      
      resListSep[[length(resListSep)+1]] <- resultsSIF
      
      rm(resultsSIF)
      
      file.remove("cplex.log")
      file.remove("results1.txt")
      file.remove("results2.txt")
      file.remove("testFile.lp")
      file.remove("testFile2.lp")
      file.remove("clone0.log")
      file.remove("clone1.log")
      
      # temp <- resultsSIF
      
    } else{
      
      delayedAssign("do.next", {next})
      
      experiments <- list(cond1=c(rownames(GMM[[1]])[i]))
      res <- tryCatch({dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=experiments)},
                      warning = function(w) {
                        force(do.next)
                      }, error = function(e) {
                        force(do.next)
                      })
      
      data.P<-dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=experiments)
      show(data.P)
      
      speciesP(data.P)
      
      res = tryCatch({build_Nw(data.On=data.P, targets.On=targets.P, bg=bg, nK = "yes")},
                     warning = function(w) {
                       force(do.next)
                     }, error = function(e) {
                       force(do.next)
                     })
      
      pknList <- buildNetwork(experiments=experiments, data.P=data.P, bg=bg, nK = "all", temp = temp)
      
      pknListTemp <- pknList
      
      show(pknList)
      
      ################################################################################
      #                         First step optimization                              #
      ################################################################################
      # Build the matrix wth the necessary data for all the species in the prior knowledge
      dataMatrix <- buildDataMatrix(dataGMM = dataGMM, pknList = pknListTemp, targets = targets, experiments = experiments)
      
      # SIF file for the prior knowledge interactions
      sif <- createSIF(pknList = pknListTemp)
      
      # Creating lists containing all our variables
      binary_x <- create_binary_variables_for_x_vector(dataMatrix = dataMatrix)
      binary_y <- create_binary_variables_for_y_vector(pknList = pknListTemp)
      binaries <- create_binaries(binaries_x = binary_x, binaries_y = binary_y)
      
      # Writing the objective function
      oF <- write_objective_function_2(dataMatrix = dataMatrix, binaries = binaries)
      
      # Writing the bounds and also all the vvariables
      bounds <- write_boundaries(binaries = binaries, pknList = pknList, M = 100)
      
      # Writing equality constraints
      eC <- write_equality_constraints(dataMatrix = dataMatrix, binaries = binaries)
      
      # Writing Constraint - 1
      c1 <- write_constraints_1(dataMatrix = dataMatrix, binaries = binaries, pknList = pknListTemp)
      
      # Writing Constraint - 2
      c2 <- write_constraints_2(dataMatrix = dataMatrix, binaries = binaries, pknList = pknListTemp)
      
      # Writing Constraint - 3
      c3 <- write_constraints_3(dataMatrix = dataMatrix, binaries = binaries, pknList = pknListTemp)
      
      # Writing Constraint - 4
      c4 <- write_constraints_4(dataMatrix = dataMatrix, binaries = binaries, pknList = pknListTemp)
      
      # Writeing Constraint - 5
      c5 <- write_constraints_5(dataMatrix = dataMatrix, binaries = binaries, pknList = pknListTemp)
      
      # Writeing Constraint - 6
      c6 <- write_self_activating_constraints(pknList = pknList, binaries = binaries, dataMatrix = dataMatrix, M = 100)
      
      # Writing Constrant - 7
      c7 <- write_time_point_constraints(binaries = binaries, tempSIF = temp, resultsSIF1 = sif)
      
      # Putting all constraints together in one file
      allC <- all_constraints(equalityConstraints = eC, constraints1 = c1, constraints2 = c2, constraints3 = c3, constraints4 = c4, constraints5 = c5, constraints6 = c7)
      
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
      temp <- unique(rbind(temp, resultsSIF1))
      write.table(resultsSIF1, file = "resultsSIF1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
      resultsSIF1 <- removeRedundantNodes(resultsSIF1 = resultsSIF1)
      resultsSIF <- unique(resultsSIF1)
      colnames(resultsSIF) <- c("Source", "Interaction", "Target")
      write.table(resultsSIF, file = paste0("resultSIF_tp", i, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
      
      resListSep[[length(resListSep)+1]] <- resultsSIF
      
      rm(resultsSIF)
      
      file.remove("cplex.log")
      file.remove("results1.txt")
      file.remove("results2.txt")
      file.remove("testFile.lp")
      file.remove("testFile2.lp")
      file.remove("clone0.log")
      file.remove("clone1.log")
      
      # temp <- unique(rbind(temp, resultsSIF))
      # temp <- resultsSIF
      
    }
    
    # write.table(temp, file = paste0("finalModelDT.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
    
  }
  
  resList[[length(resList)+1]] <- resListSep
  
  iterSum <- iterSum + 1
  print(paste0("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@____________",iterSum))
  
}

save(resList, file = "resList.RData")

end.time<-Sys.time()
end.time-start.time






