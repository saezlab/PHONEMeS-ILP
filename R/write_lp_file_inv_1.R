#'\code{write_lp_file_inv_1}
#'
#'Enio Gjerga
#'

write_lp_file_inv_1 <- function(dataGMM, pknList, targets, experiments, penFac){
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
  oF <- write_objective_function(dataMatrix = dataMatrix, binaries = binaries, sizePen = penFac)
  
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
  # write(paste0("set mip tolerances mipgap ", 0), data2, append = TRUE)
  # write(paste0("set mip pool relgap ", 0), data2, append = TRUE)
  write("optimize", data2, append = TRUE)
  write("write results1.txt sol all", data2, append = TRUE)
  write("quit", data2, append = TRUE)
  
}