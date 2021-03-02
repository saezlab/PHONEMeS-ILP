\name{runPHONEMeS}
\alias{runPHONEMeS}
\title{
Function to run PHONEMeS-ILP
}
\description{
This function takes the PHONEMeS inputs (data object, background network, drug perturbation targets, experimental design and solver parameters) and returns the network solution
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
  \item{nSolutions}{
number of solution to be returned by PHONEMeS (default: 100)
}
  \item{mipgap}{
CPLEX/Cbc parameter: the absolute tolerance on the gap between the best integer objective and the objective of the best node remaining. When this difference falls below the value of this parameter, the linear integer optimization is stopped. Default set to 0
}
  \item{relgap}{
CPLEX/Cbc parameter: Allowed relative gap of accepted solution comparing within the pool of accepted solution (fraction; default: 0)
}
  \item{replace}{
CPLEX parameter: Replacement strategy of solutions in the pool (0,1,2 - default: 2 = most diversified solutions)
}
  \item{populate}{
CPLEX parameter: Allowed number of solutions to be generated (default: 5000)
}
  \item{intensity}{
CPLEX parameter: Intensity of solution searching (0,1,2,3,4 - default: 4)
}
  \item{timelimit}{
CPLEX/Cbc parameter: Time limit of CPLEX optimisation (in seconds - default: 3600)
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
