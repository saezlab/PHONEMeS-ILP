# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("buildDataMatrix.R")
source("ilpFunctions.R")

load(file='allD_noCSK_filt.RData')
load(file='data_Hsu_PHONEMeS.RData')

#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)
#Choose the drug targets
targets <- list(cond1=c("MTOR_HUMAN"))
targets.P<-targets
#Choose the drug treatments matching to the drug targets
#and match to what is present in the background network
experiments <- list(cond1=c("tor - ins"))
data.P<-dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=experiments)
show(data.P)

speciesP(data.P)

#Create the PKN list that will be used for optimisation
pknList<-buildNw(data.On=data.P, targets.On=targets.P, bg=bg,nK="no")
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

system(paste0(getwd(), "/cplex -f cplexCommand.txt"))


# Read the results from the CPLEX and do the necessary processing of the model
library(XML)
resultsSIF1 <- readOutSIF(cplexSolutionFileName = "results1.txt", binaries = binaries)
colnames(resultsSIF1) <- c("Source", "Interaction", "Target")
write.table(resultsSIF1, file = "resultsSIF.txt", quote = FALSE, row.names = FALSE, sep = "\t")

file.remove("cplex.log")
file.remove("results1.txt")
file.remove("testFile.lp")

resultsSIF <- resultsSIF1
