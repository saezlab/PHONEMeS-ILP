\name{runLimma}
\alias{runLimma}
\title{
Function to perform differential analysis with Limma
}
\description{
This the phospho-proteomics measurement data and the experimental design of conditions for each sample and performs differential analysis between groups of samples.
}
\arguments{
  \item{measurements}{
phosphoproteomics measurements dataset
}
  \item{targets}{
matrix indicationg the conditions assigned to each sample	
}
  \item{comparisons}{
list object indicating which groups of samples to compare
}
  \item{pool}{
pool
}
  \item{regress_out}{
regress_out
}
}
\value{
This function returns the differential analysis results.
}
\author{
A. Dugourd
}
