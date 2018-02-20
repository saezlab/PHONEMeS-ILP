runPHONEMeSUD <- function(dataInput = dataInput, bg = bg, targets.P = targets.P, data.P = data.P, experiments = experiments){
  
  conditions <- experiments
  
  resultsSIF <- matrix(data = , nrow = 1, ncol = 3)
  
  #Create the PKN list that will be used for optimisation
  pknList<-buildNw(data.On=data.P, targets.On=targets.P, bg=bg,nK="all")
  
  show(pknList)
  
  targets <- unlist(targets.P)
  
  # Build the matrix wth the necessary data for all the species in the prior knowledge
  dataMatrix <- buildDataMatrix(dataGMM = dataGMM, pknList = pknList, targets = targets, experiments = experiments)
  
  # SIF file for the prior knowledge interactions
  sif <- createSIF(pknList = pknList)
  
  # Creating lists containing all our variables
  binary_x <- create_binary_variables_for_x_vector(dataMatrix = dataMatrix)
  binary_y <- create_binary_variables_for_y_vector(pknList = pknList)
  binaries <- create_binaries(binaries_x = binary_x, binaries_y = binary_y)
  
  # Writing the objective function
  oF <- write_objective_function(dataMatrix = dataMatrix, binaries = binaries, sizePen = TRUE, sizePenType = "edge", lambda = 0.001)
  
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
  c6 <- write_self_activating_constraints(sif = createSIF(pknList = pknList), binaries = binaries, dataMatrix = dataMatrix, M = 100)
  
  # Putting all constraints together in one file
  allC <- all_constraints(equalityConstraints = eC, constraints1 = c1, constraints2 = c2, constraints3 = c3, constraints4 = c4, constraints5 = c5, constraints6 = c6)
  
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
  resultsSIF1 <- readOutSIF(cplexSolutionFileName = "results1.txt", binaries = binaries, pknList = pknList, dataMatrix = dataMatrix)
  colnames(resultsSIF1) <- c("Source", "Interaction", "Target")
  write.table(unique(resultsSIF1), file = "resultsSIF-Downstream.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  resultsSIFDownstream <- unique(resultsSIF1)
  resultsSIF <- rbind(resultsSIF, unique(resultsSIF1))[-1, ]
  
  file.remove("cplex.log")
  for (i in 0:20)
  {
    file.remove(paste("clone",paste(i,".log", sep = ""), sep = ""))
  }
  file.remove("results1.txt")
  file.remove("testFile.lp")
  
  #############################
  ##  Downside Modelling    ##
  #############################
  
  data.P <- dataBycond(dataGMM, bg, scaled = TRUE, rowBycond = conditions)
  
  show(data.P)
  
  speciesP(data.P)
  #Create the PKN list that will be used for optimisation
  pknList<-build_Nw_Inv(data.On=data.P, targets.On=targets.P, bg=bg,nK="all")
  
  show(pknList)
  
  targets <- unlist(targets.P)
  
  # Build the matrix wth the necessary data for all the species in the prior knowledge
  dataMatrix <- buildDataMatrix(dataGMM = dataGMM, pknList = pknList, targets = targets, experiments = experiments)
  
  # SIF file for the prior knowledge interactions
  sif <- createSIF(pknList = pknList)
  
  # Creating lists containing all our variables
  binary_x <- create_binary_variables_for_x_vector(dataMatrix = dataMatrix)
  binary_y <- create_binary_variables_for_y_vector(pknList = pknList)
  binaries <- create_binaries(binaries_x = binary_x, binaries_y = binary_y)
  
  # Writing the objective function
  oF <- write_objective_function(dataMatrix = dataMatrix, binaries = binaries, sizePen = TRUE, sizePenType = "edge", lambda = 0.001)
  
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
  c6 <- write_self_activating_constraints(sif = createSIF(pknList = pknList), binaries = binaries, dataMatrix = dataMatrix, M = 100)
  
  # Putting all constraints together in one file
  allC <- all_constraints(equalityConstraints = eC, constraints1 = c1, constraints2 = c2, constraints3 = c3, constraints4 = c4, constraints5 = c5, constraints6 = c6)
  
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
  resultsSIF1 <- readOutSIF(cplexSolutionFileName = "results1.txt", binaries = binaries, pknList = pknList, dataMatrix = dataMatrix)
  colnames(resultsSIF1) <- c("Target", "Interaction", "Source")
  resultsSIF1 <- resultsSIF1[,c(3,1,2)]
  
  write.table(unique(resultsSIF1), file = "resultsSIF-Upstream.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  resultsSIFUpstream <- unique(resultsSIF1)
  resultsSIF1 <- resultsSIF1[,c(1,3,2)]
  resultsSIF <- unique(rbind(resultsSIF, unique(resultsSIF1)))
  write.table(unique(resultsSIF)[-1, ], file = "resultsSIF-UD.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  file.remove("cplex.log")
  for (i in 0:20)
  {
    file.remove(paste("clone",paste(i,".log", sep = ""), sep = ""))
  }
  file.remove("results1.txt")
  file.remove("testFile.lp")
  
  species <- unique(c(resultsSIF[, 1], resultsSIF[, 3]))
  measuredSites <- species[which(species%in%dataInput$IDmap$S.cc)]
  targets <- gsub(pattern = "TG:", replacement = "", x = colnames(dataMatrix$dataMatrix)[dataMatrix$tgID])
  
  nodesAttributes <- matrix(data = , nrow = length(species), ncol = 2)
  colnames(nodesAttributes) <- c("species", "nodesP")
  
  nodesAttributes[, 1] <- species
  nodesAttributes[which(nodesAttributes[, 1]%in%targets), 2] <- "D"
  nodesAttributes[which(nodesAttributes[, 1]%in%measuredSites), 2] <- "P"
  
  write.table(nodesAttributes, file = "nodesAttributes-UD.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  result <- list()
  result[[1]] <- resultsSIFDownstream
  result[[2]] <- resultsSIFUpstream
  result[[3]] <- unique(resultsSIF)[-1, ]
  result[[4]] <- nodesAttributes
  
  names(result) <- c("Downstream", "Upstream", "Combined", "nodesAttributes")
  
  return(result)
  
  
}