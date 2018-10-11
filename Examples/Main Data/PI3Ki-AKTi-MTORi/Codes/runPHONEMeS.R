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

#' Run PHONEMeS ILP
#' 
#' @param Arguments 
#' @param targets.P
#' @param conditions
#' @param dataGMM
#' @param experiments
#' @param bg
#' @param nK 
#' @param solver Solver to use for solving the ILP.
#
#' @return SIF like data.frame with the output network.
runPHONEMeS <- function(targets.P, conditions, dataGMM, experiments, bg, nK="all", solver="cplex"){
  
  conditions <- conditions[experiments]
  valid_solver_list <- c("cplex", "cbc")
  if (!(solver %in% valid_solver_list)){
    stop(paste0("Select a valid solver option (", paste(valid_solver_list, collapse=", "), ")"))
  }
  
  data.P <- dataBycond(dataGMM, bg, scaled=TRUE, rowBycond=conditions)
  show(data.P)
  
  speciesP(data.P)
  
  targets.P <- targets.P[experiments]
  
  pknList<-build_Nw(data.On=data.P, targets.On=targets.P, bg=bg, nK=nK)
  
  show(pknList)
  
  TG <- unique(unlist(targets.P))
  
  if(length(TG)==1){
    
    targets <- targets.P
    
    write_lp_file(dataGMM = dataGMM, pknList = pknList, targets = targets, experiments = conditions[experiments])
    
    if (solver=="cplex"){
      resultsSIF1 <- solve_with_cplex()
    } else if (solver=="cbc"){
      resultsSIF1 <- solve_with_cbc()
    } else {
      stop("Select a valid solver option ('cplex', 'cbc')")
    }
    
    resultsSIF <- resultsSIF1
    
  } else {
    
    for(ii in 1:length(TG)){
      
      idxT <- c()
      for(i in 1:length(targets.P)){
        if(TG[ii]%in%targets.P[[i]]){
          idxT <- c(idxT, i)
        }
      }
      
      write_lp_file(dataGMM = dataGMM, pknList = pknList, targets = targets.P[idxT], experiments = conditions[idxT])
      
      if (solver=="cplex"){
        resultsSIF1 <- solve_with_cplex()
      } else if (solver=="cbc"){
        resultsSIF1 <- solve_with_cbc()
      } else {
        stop("Select a valid solver option ('cplex', 'cbc')")
      }
      
      if(ii==1){
        resultsSIF <- resultsSIF1
      } else {
        resultsSIF <- unique(rbind(resultsSIF, resultsSIF1))
      }
      
    }
    
  }
  
  if(length(which(duplicated(resultsSIF[, c(1, 3)]))) > 0){
    
    returnSIF <- resultsSIF[-which(duplicated(resultsSIF[, c(1, 3)])), ]
    
    for(ii in nrow(returnSIF)){
      
      idx1 <- which(resultsSIF[, 1]==returnSIF[ii, 1])
      idx2 <- which(resultsSIF[, 3]==returnSIF[ii, 3])
      
      idx <- intersect(x = idx1, y = idx2)
      
      returnSIF[ii, 2] <- as.character(mean(as.numeric(resultsSIF[idx, 2])))
      
    }
    
  } else {
    
    returnSIF <- resultsSIF
    
  }
  
  colnames(returnSIF) <- c("Source", "Weight", "Target")
  
  return(returnSIF)
  
}


write_lp_file <- function(dataGMM, pknList, targets, experiments){
  # This function writes the optimization problem to be solved in a *.lp file
  
  # Build the matrix wth the necessary data for all the species in the prior knowledge
  dataMatrix <- buildDataMatrix(dataGMM = dataGMM, pknList = pknList, targets = targets, experiments = experiments)
  
  # SIF file for the prior knowledge interactions
  sif <- createSIF(pknList = pknList)
  
  # Creating lists containing all our variables
  binary_x <- create_binary_variables_for_x_vector(dataMatrix = dataMatrix)
  binary_y <- create_binary_variables_for_y_vector(pknList = pknList)
  binaries <- create_binaries(binaries_x = binary_x, binaries_y = binary_y)
  
  # save mapping information
  saveRDS(binaries, file="tmp_binaries.rds")
  
  # Writing the objective function
  oF <- write_objective_function(dataMatrix = dataMatrix, binaries = binaries, sizePen = TRUE, penMode = "edge")
  
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
  
  # Putting all constraints together in one file
  allC <- all_constraints(equalityConstraints = eC, constraints1 = c1, constraints2 = c2, constraints3 = c3, constraints4 = c4, constraints5 = c5, constraints6 = c6)
  
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
}


solve_with_cplex <- function(){
  system(paste0(getwd(), "/cplex -f cplexCommand.txt"))
  
  # load mapping information
  binaries <- readRDS("tmp_binaries.rds")
  
  # Read the results from the CPLEX and do the necessary processing of the model
  library(XML)
  resultsSIF1 <- readOutSIFAll(cplexSolutionFileName = "results1.txt", binaries = binaries)
  if(length(which(duplicated(resultsSIF1[, c(1, 3)]))) > 0){
    
    returnSIF <- resultsSIF1[-which(duplicated(resultsSIF1[, c(1, 3)])), ]
    
    for(ii in nrow(returnSIF)){
      
      idx1 <- which(resultsSIF1[, 1]==returnSIF[ii, 1])
      idx2 <- which(resultsSIF1[, 3]==returnSIF[ii, 3])
      
      idx <- intersect(x = idx1, y = idx2)
      
      returnSIF[ii, 2] <- as.character(mean(as.numeric(resultsSIF1[idx, 2])))
      
    }
    
  } else {
    
    returnSIF <- resultsSIF1
    
  }
  
  colnames(returnSIF) <- c("Source", "Weight", "Target")
  
  resultsSIF1 <- returnSIF
  # change format to data.frame
  resultsSIF1 <- data.frame(resultsSIF1, stringsAsFactors=FALSE)
  resultsSIF1 <- resultsSIF1 %>% mutate(Weight=as.numeric(Weight))
  
  # clean-up temporary files
  file.remove("cplex.log")
  file.remove("results1.txt")
  file.remove("testFile.lp")
  file.remove("tmp_binaries.rds")
  
  return(resultsSIF1)
}


solve_with_cbc <- function(){
  cbc_command <- "cbc testFile.lp solve printi csv solu results_cbc.txt"
  system(cbc_command)
  
  # retrieve solution
  readCbcSolution <- function(file, binaries){
    library("dplyr")
    library("tidyr")
    
    # read cbc solution file
    cbc_table <- read.csv("results_cbc.txt")
    
    # load mapping information
    binaries <- readRDS("tmp_binaries.rds")
    
    # find mapping of variables from the MILP to the model
    mapping_table <- data.frame(binaries)
    colnames(mapping_table) <- c("milp", "math", "description")
    # keep only reaction variables
    mapping_table <- mapping_table %>% filter(grepl("reaction", description))
    # merge with solution
    cbc_table <- merge(cbc_table, mapping_table, by.x="name", by.y="milp")
    # keep only active reactions
    cbc_table <- cbc_table %>% mutate(solution=round(as.numeric(solution))) %>% 
      filter(solution==1)
    # create SIF table
    sif <- cbc_table %>% select(description) %>% 
      mutate(description = gsub("reaction ", "", description)) %>% 
      separate(description, into=c("Source", "Target"), sep="=") %>% 
      mutate(Weight = 1) %>% 
      select(Source, Weight, Target)
    
    return(sif)
  }
  
  resultsSIF1 <- readCbcSolution("results_cbc.txt", binaries)
  
  # clean-up temporary files
  file.remove("results_cbc.txt")
  file.remove("testFile.lp")
  file.remove("tmp_binaries.rds")
  
  return(resultsSIF1)
}
