\name{runPHONEMeS_mult}
\alias{runPHONEMeS_mult}
\title{
Function to run the PHONEMeS-ILP method for time-resolved phosphoproteomics data. 
}
\description{
This function takes the PHONEMeS inputs of the time-point variant (data object, background network, drug perturbation targets, experimental design dependent on the time-points and solver parameters) and returns the time-resolved network solution
}
\arguments{
  \item{targets.P}{
a list of drug targets
}
  \item{conditions}{
the experimental conditions from the data	
}
  \item{inputObj}{
data inputs of PHONEMeS
}
  \item{experiments}{
which experimental conditions to consider for the analysis
}
  \item{bg}{
the background network
}
  \item{nK}{
a single string ("all","no", "drugs2data", "data") specifying whether the predicted edges should be included at all levels, none, only in the network that connects drug targets to data kinases or only between data sites and their kinases
}
  \item{solver}{
solver to use (either: cplex or cbc - cplex by default)
}
  \item{solverPath}{
path to solver (default: "~/Documents/cplex")
}
  \item{nIter}{
number of iterations to run/solution to be returned by PHONEMeS (default: 100)
}
  \item{penFac}{
size penalty factor (default: 0.00001)
}
}
\value{
This function returns the combined network solutions of the time-resolved variant of PHONEMeS. It contains the Source, the Target, the Weight and the Time-Point. The Weight tells how frequently an interaction is present across the different solutioins. The Time-Point tells after which time-point an interaction appears in the solution.
}
\author{
E. Gjerga
}
