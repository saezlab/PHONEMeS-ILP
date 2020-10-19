\name{buildInputs}
\alias{buildInputs}
\title{
Function to generate the data inputs of PHONEMeS
}
\description{
Building the PHONEMeS inputs from the differential analysis results obtained from the Limma analysis
}
\arguments{
  \item{tableTopList}{
a list of table-top contiaing the differential analysis results obtained from Limma. Tables for each comparison/experiment should be combined in one list.
}
  \item{fcThresh}{
a numberc value or a vector of values (for each experimental condition) indicating the threshold of significance (as a log2) (sites with regulation levels higher than the set significance are considered as perturbed. default set to NULL)	
}
  \item{pThresh}{
a p-value threshold parameter (or vector if different threshold should be applied for different comparisons) to indicate which sites/measurements are considered as significant
}
  \item{idxID}{
the index indicating which column of the table-top contains the sites identifiers (default: 1)
}
  \item{idxFC}{
the index indicating which column of the table-top contains the FC values (default: 2)
}
  \item{idxPval}{
the index indicating which column of the table-top contains the p-values (default: 6, the adjusted p-values column of the table-tops returned by Limma)
}
  \item{mappingTable}{
if user wishes to plot on the solution identifiers different than what is on the measurements table then he should provide a mapping table (default: NULL)
}
  \item{namesConditions}{
a vector of experimental conditions to name each experiment in hte output
}
  \item{directionality}{
a string character indicating whether the directionality of the measurement should be considered for a site to be consiidered as perturbed. It can take three values: up - measurements with log2FC which are positive and higher than fcThresh are considered as perturbed; down - measurements with log2FC which are negative and lower than fcThresh are considered as perturbed; no - measurement whose absolute log2FC value is higher than fcThresh are considered as perturbed. Default: no.
}
}
\value{
This function returns the data inputs used for the PHONEMeS analysis.
}
\author{
E. Gjerga
}
