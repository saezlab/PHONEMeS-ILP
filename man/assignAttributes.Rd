\name{assignAttributes}
\alias{assignAttributes}
\title{
Assigning node attributes
}
\description{
This function assigns attributes to each node in the PHONEMeS solution (D - whether a node is a drug target; P - whther a node is a measured site; while the other nodes which have no attribute assigned are the ones that are inferred by PHONEMeS). This function is used anly to assist the visualization of the PHONEMeS solutions in CytoScape.
}
\arguments{
  \item{sif}{
the PHONEMeS solution as a matrix of interactions with three columns (Target, Source and Weight)
}
  \item{inputObj}{
the PHONEMeS data inputs that were used for the analysis	
}
  \item{targets}{
the list of drug perturbation targets
}
  \item{writeAttr}{
a verbose (TRUE/FALSE) indicating whether to write the attributes as a txt file
}
  \item{fileName}{
the name of the file and the path where to save the txt files of the attributes
}
}
\value{
This function returns a matrix with two columns where in the first column we have the node ID's and on the second column we have the attributes assigned (D, P or none)
}
\author{
E. Gjerga
}
