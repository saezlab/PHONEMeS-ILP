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
runPHONEMeS_UD <- function(targets.P, conditions, dataGMM, experiments, bg, nK="all", solver="cplex", mipgap=0, relgap=0, timelimit=300){
  
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
  
  write_lp_file_1(dataGMM = dataGMM, pknList = pknList, targets = targets.P, experiments = conditions, mipgap=mipgap, relgap=relgap, timelimit=timelimit)
  
  if (solver=="cplex"){
    resultsSIF1 <- solve_with_cplex()
  } else if (solver=="cbc"){
    resultsSIF1 <- solve_with_cbc()
  } else {
    stop("Select a valid solver option ('cplex', 'cbc')")
  }
  
  # write.table(resultsSIF, file = "resultsSIF.txt", quote = FALSE, row.names = FALSE, sep = "\t")
  colnames(resultsSIF1) <- c("Source", "Interaction", "Target")
  resultsSIF1[, 2] <- "1"
  
  #####
  # Step - 2
  #Create the PKN list that will be used for optimisation
  pknList <- build_Nw_Inv(data.On=data.P, targets.On=targets.P, bg=bg,nK=nK)
  
  if(!is.null(pknList)){
    show(pknList)
    
    TG <- unlist(targets.P)
    
    write_lp_file_2(dataGMM = dataGMM, pknList = pknList, targets = targets.P, experiments = conditions, mipgap=mipgap, relgap=relgap, timelimit=timelimit)
    
    if (solver=="cplex"){
      resultsSIF2 <- solve_with_cplex()
    } else if (solver=="cbc"){
      resultsSIF2 <- solve_with_cbc()
    } else {
      stop("Select a valid solver option ('cplex', 'cbc')")
    }
    
    # write.table(resultsSIF, file = "resultsSIF.txt", quote = FALSE, row.names = FALSE, sep = "\t")
    colnames(resultsSIF2) <- c("Source", "Interaction", "Target")
    resultsSIF2[, 2] <- "1"
    
    resultSIF <- resultsSIF1
    for(ii in 1:nrow(resultsSIF2)){
      
      ss <- resultsSIF2[ii, 1]
      tt <- resultsSIF2[ii, 3]
      idx1 <- which(resultSIF[, 1]==ss)
      idx2 <- which(resultSIF[, 3]==tt)
      idx <- intersect(x = idx1, y = idx2)
      
      if(length(idx)>0){
        resultSIF[idx, 2] <- mean(resultSIF[idx, 2], resultsSIF2[ii, 2])
      } else {
        resultSIF <- rbind(resultSIF, resultsSIF2[ii, ])
      }
      
    }
    
    resList <- list()
    resList[[1]] <- resultsSIF1
    resList[[2]] <- resultsSIF2
    resList[[3]] <- resultSIF
    
    names(resList) <- c("Downside", "Upside", "Combined")
    
    return(resList)
    
  } else {
    
    resList <- list()
    resList[[1]] <- resultsSIF1
    resList[[2]] <- NULL
    resultSIF <- resultsSIF1
    resList[[3]] <- resultSIF
    
    return(resList)
    
  }
  
}


write_lp_file_1 <- function(dataGMM, pknList, targets, experiments, mipgap, relgap, timelimit){
  # This function writes the optimization problem to be solved in a *.lp file
  
  # Build the matrix wth the necessary data for all the species in the prior knowledge
  dataMatrix <- buildDataMatrix(dataGMM = dataGMM, pknList = pknList, targets = targets, experiments = experiments)
  
  # SIF file for the prior knowledge interactions
  sif <- createSIF(pknList = pknList)
  
  # Creating lists containing all our variables
  binary_x <- create_binary_variables_for_x_vector(dataMatrix = dataMatrix)
  binary_y <- create_binary_variables_for_y_vector(pknList = pknList)
  binary_z <- create_binary_variables_for_z_vector(pknList = pknList, dataMatrix = dataMatrix)
  binary_orin_orout <- create_orin_orout_variables(dataMatrix = dataMatrix, targets = targets, pknList = pknList)
  binaries <- create_binaries(binaries_x = binary_x, binaries_z = binary_z, binaries_in_out = binary_orin_orout, binaries_y = binary_y)
  
  # save mapping information
  saveRDS(binaries, file="tmp_binaries.rds")
  
  # Writing the objective function
  oF <- write_objective_function(dataMatrix = dataMatrix, binaries = binaries)
  
  # Writing the bounds and also all the vvariables
  bounds <- write_boundaries(binaries = binaries, pknList = pknList, M = 100, dataMatrix = dataMatrix)
  
  # Writing equality constraints
  eC <- write_equality_constraints(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
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
  
  # Writeing Constraint - 7
  c7 <- write_in_out_constraints(binaries = binaries, targets = targets, dataMatrix = dataMatrix, pknList = pknList)
  
  # Putting all constraints together in one file
  allC <- all_constraints(equalityConstraints = eC, constraints1 = c1, constraints2 = c2, constraints3 = c3, constraints4 = c4, constraints5 = c5, constraints6 = c(c6, c7))
  
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
  
  # write cplexCommand file
  # data2 = "cplexCommand.txt"
  # write("read testFile.lp", data2)
  # write(paste0("set mip tolerances mipgap ", mipgap), data2, append = TRUE)
  # write(paste0("set mip pool relgap ", relgap), data2, append = TRUE)
  # write(paste0("set mip pool replace ", replace), data2, append = TRUE)
  # write(paste0("set mip limits populate ", populate), data2, append = TRUE)
  # write(paste0("set mip pool capacity ", nSolutions), data2, append = TRUE)
  # write(paste0("set mip pool intensity ", intensity), data2, append = TRUE)
  # write(paste0("set timelimit ", timelimit), data2, append = TRUE)
  # write("populate", data2, append = TRUE)
  # write("write results1.txt sol all", data2, append = TRUE)
  # write("quit", data2, append = TRUE)
  
  data2 = "cplexCommand.txt"
  write("read testFile.lp", data2)
  write(paste0("set mip tolerances mipgap ", 0), data2, append = TRUE)
  write(paste0("set mip pool relgap ", 0), data2, append = TRUE)
  write("optimize", data2, append = TRUE)
  write("write results1.txt sol all", data2, append = TRUE)
  write("quit", data2, append = TRUE)
  
}

write_lp_file_2 <- function(dataGMM, pknList, targets, experiments, mipgap, relgap, timelimit){
  # This function writes the optimization problem to be solved in a *.lp file
  
  # Build the matrix wth the necessary data for all the species in the prior knowledge
  dataMatrix <- buildDataMatrix(dataGMM = dataGMM, pknList = pknList, targets = targets, experiments = experiments)
  
  # SIF file for the prior knowledge interactions
  sif <- createSIF(pknList = pknList)
  
  # Creating lists containing all our variables
  binary_x <- create_binary_variables_for_x_vector(dataMatrix = dataMatrix)
  binary_y <- create_binary_variables_for_y_vector(pknList = pknList)
  binary_z <- create_binary_variables_for_z_vector(pknList = pknList, dataMatrix = dataMatrix)
  binary_orin_orout <- create_orin_orout_variables(dataMatrix = dataMatrix, targets = targets, pknList = pknList)
  binaries <- create_binaries(binaries_x = binary_x, binaries_z = binary_z, binaries_in_out = binary_orin_orout, binaries_y = binary_y)
  
  # save mapping information
  saveRDS(binaries, file="tmp_binaries.rds")
  
  # Writing the objective function
  oF <- write_objective_function(dataMatrix = dataMatrix, binaries = binaries)
  
  # Writing the bounds and also all the vvariables
  bounds <- write_boundaries(binaries = binaries, pknList = pknList, M = 100, dataMatrix = dataMatrix)
  
  # Writing equality constraints
  eC <- write_equality_constraints(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Writing Constraint - 1
  c1 <- write_constraints_1(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Writing Constraint - 2
  c2 <- write_constraints_2_ud(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Writing Constraint - 3
  c3 <- write_constraints_3(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Writing Constraint - 4
  c4 <- write_constraints_4(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Writeing Constraint - 5
  c5 <- write_constraints_5_ud(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList)
  
  # Writeing Constraint - 6
  c6 <- write_self_activating_constraints(pknList = pknList, binaries = binaries, dataMatrix = dataMatrix, M = 100)
  
  # Writeing Constraint - 7
  c7 <- write_in_out_constraints(binaries = binaries, targets = targets, dataMatrix = dataMatrix, pknList = pknList)
  
  # Putting all constraints together in one file
  allC <- all_constraints(equalityConstraints = eC, constraints1 = c1, constraints2 = c2, constraints3 = c3, constraints4 = c4, constraints5 = c5, constraints6 = c(c6, c7))
  
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
  
  # write cplexCommand file
  # data2 = "cplexCommand.txt"
  # write("read testFile.lp", data2)
  # write(paste0("set mip tolerances mipgap ", mipgap), data2, append = TRUE)
  # write(paste0("set mip pool relgap ", relgap), data2, append = TRUE)
  # write(paste0("set mip pool replace ", replace), data2, append = TRUE)
  # write(paste0("set mip limits populate ", populate), data2, append = TRUE)
  # write(paste0("set mip pool capacity ", nSolutions), data2, append = TRUE)
  # write(paste0("set mip pool intensity ", intensity), data2, append = TRUE)
  # write(paste0("set timelimit ", timelimit), data2, append = TRUE)
  # write("populate", data2, append = TRUE)
  # write("write results1.txt sol all", data2, append = TRUE)
  # write("quit", data2, append = TRUE)
  
  data2 = "cplexCommand.txt"
  write("read testFile.lp", data2)
  write(paste0("set mip tolerances mipgap ", 0), data2, append = TRUE)
  write(paste0("set mip pool relgap ", 0), data2, append = TRUE)
  write("optimize", data2, append = TRUE)
  write("write results1.txt sol all", data2, append = TRUE)
  write("quit", data2, append = TRUE)
  
}


solve_with_cplex <- function(){
  system(paste0(getwd(), "/cplex -f cplexCommand.txt"))
  
  # load mapping information
  binaries <- readRDS("tmp_binaries.rds")
  
  # Read the results from the CPLEX and do the necessary processing of the model
  library(XML)
  resultsSIF1 <- readOutSIFAll(cplexSolutionFileName = "results1.txt", binaries = binaries)
  colnames(resultsSIF1) <- c("Source", "Interaction", "Target")
  # write.table(resultsSIF1, file = "res1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
  
  # change format to data.frame
  resultsSIF1 <- data.frame(resultsSIF1, stringsAsFactors=FALSE)
  resultsSIF1 <- resultsSIF1 %>% mutate(Interaction=as.numeric(Interaction))
  
  # clean-up temporary files
  file.remove("cplex.log")
  for (i in 0:20)
  {
    file.remove(paste("clone",paste(i,".log", sep = ""), sep = ""))
  }
  file.remove("results1.txt")
  file.remove("testFile.lp")
  file.remove("tmp_binaries.rds")
  file.remove("cplexCommand.txt")
  
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
      mutate(Interaction = 1) %>% 
      select(Source, Interaction, Target)
    
    return(sif)
  }
  
  resultsSIF1 <- readCbcSolution("results_cbc.txt", binaries)
  
  # clean-up temporary files
  file.remove("results_cbc.txt")
  file.remove("testFile.lp")
  file.remove("tmp_binaries.rds")
  
  return(resultsSIF1)
}