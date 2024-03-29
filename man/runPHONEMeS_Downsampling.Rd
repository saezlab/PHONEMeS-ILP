\name{runPHONEMeS_Downsampling}
\alias{runPHONEMeS_Downsampling}
\title{
Function to run PHONEMeS-ILP. Similar to runPHONEMeS, however the multiple solutions in this case are generated by running the PHONEMeS analysis multiple times with random downsampling of the data and then the solutions are combined in one consensus network.
}
\description{
This function takes the PHONEMeS inputs and returns the network solution
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
path to solver (default: "/usr/bin/cplex")
}
  \item{nIter}{
number of iterations to run/solution to be returned by PHONEMeS (default: 100)
}
  \item{penFac}{
size penalty factor (default: 0.00001)
}
}
\value{
This function returns the combined network solutions of PHONEMeS. It contains the Source, the Target and the Weight of each interaction. The Weight tells how frequently an interaction is present across the different solutioins.
}
\author{
E. Gjerga
}
