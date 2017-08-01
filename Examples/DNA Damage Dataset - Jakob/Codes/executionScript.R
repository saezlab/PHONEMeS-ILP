start.time <- Sys.time()
#load the packages
library(BioNet)
library(igraph)
library(PHONEMeS)

load("data4cluster_1.RData")
targets <- targets.P
experiments <- list(cond1=c("ATM - KU60019"))

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
oF2 <- write_objective_function_step_2(binaries = binaries2, dataMatrix = dataMatrix2)
bounds2 <- write_boundaries_step_2(binaries = binaries2)
c1 <- write_fixed_nodes_constraints(dataMatrix = dataMatrix2, binaries = binaries2)
c2 <- write_one_edge_constraint(sif = resultsSIF1, dataMatrix = dataMatrix2, binaries = binaries2)
c3 <- write_intermediate_node_constraints_in(sif = resultsSIF1, dataMatrix = dataMatrix2, binaries = binaries2)
c4 <- write_intermediate_node_constraints_out(sif = resultsSIF1, dataMatrix = dataMatrix2, binaries = binaries2)
c5 <- write_reaction_constraints(dataMatrix = dataMatrix2, sif = resultsSIF1, binaries = binaries2)
allC <- all_constraints_2(constraints1 = c1, constraints2 = c2, constraints3 = c3, constraints4 = c4, constraints5 = c5)

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

resultsSIF2 <- readOutSIF(cplexSolutionFileName = "results2.txt", binaries = binaries2)
resultsSIF <- reducedModel(sif = resultsSIF2, dataMatrix = dataMatrix)
write.table(resultsSIF, file = "resultsSIF.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

