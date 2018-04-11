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
load(file='dataObjects_PHONEMeS.RData')

#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM.res.noFC, IDmap=GMM.res.ID, resFC=GMM.res)

conditions <- list(c("AKT1 - Control", "AKT2 - Control"), c("CAMK1 - Control", "CAMK2 - Control"),
                    c("EGFR1 - Control", "EGFR2 - Control"), c("ERK1 - Control", "ERK2 - Control"),
                    c("MEK1 - Control", "MEK2 - Control"), c("MTOR1 - Control", "MTOR2 - Control"),
                    c("P70S6K1 - Control", "P70S6K2 - Control"), c("PI3K1 - Control", "PI3K2 - Control"),
                    c("PKC1 - Control", "PKC2 - Control"), c("ROCK1 - Control", "ROCK2 - Control"))

names(conditions) <- c("AKT1_HUMAN", "KCC2D_HUMAN", "EGFR_HUMAN", "MK01_HUMAN",
                        "MP2K1_HUMAN", "MTOR_HUMAN",  "KS6B1_HUMAN", "PK3CA_HUMAN",
                        "KPCA_HUMAN", "ROCK1_HUMAN")

targets.P<-list(cond1=c(), cond2=c(), cond3=c(), cond4=c(),
                cond5=c(), cond6=c("MTOR_HUMAN"), cond7=c(), cond8=c(),
                cond9=c(), cond10=c())

targets <- targets.P

experiments <- conditions[c(6)]
data.P<-dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=experiments)
show(data.P)

speciesP(data.P)

pknList<-buildNw(data.On=data.P, targets.On=targets.P, bg=bg, nK = "no")

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
