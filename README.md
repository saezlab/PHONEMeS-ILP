# PHONEMeS-ILP
ILP implementation of PHONEMeS - Enio GJERGA

**PHONEMeS** (**PHO**sphorylation **NE**tworks for **M**ass **S**pectrometry) is a method to model signalling networks based on untargeted phosphoproteomics mass spectrometry data and kinase/phosphatase-substrate interactions. 
Please see [Terfve et al.](http://www.nature.com/articles/ncomms9033) for an explanation of the methodolgy.

This repository contains the scripts for the ILP (Integer Linear Programming) implementation of the [PHONEMeS R package](https://github.com/saezlab/PHONEMeS/tree/master/Package) and accompanying scripts that implement the method. ILP is a mathematical optimisation algorithm in which the objective function and constraints are linear and the variables are integers.

### License

Distributed under the GNU GPLv3 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/PHONEMeS-ILP/blob/master/LICENSE).

### Installation

Before using the method, please install the current R package for [PHONEMeS](https://github.com/saezlab/PHONEMeS). For installation, download the tar file of the package and type in R:

```R
# Install PHONEMeS from Github using devtools
# install.packages('devtools') # in case devtools hasn't been installed
library(devtools)
install_github('saezlab/PHONEMeS-ILP')
# or download the source file from GitHub and install from source
install.packages('path_to_extracted_CARNIVAL_directory', repos = NULL, type="source")
```

## Running PHONEMeS

The PHONEMeS library can be initialized by:

```R
library(PHONEMeS)
```

A list of tutorials/examples of the PHONEMeS package, can be found [here](https://github.com/saezlab/PHONEMeS-ILP-Examples)

### Prerequisites

PHONEMeS requires the interactive version of IBM Cplex or CBC-COIN solver as the network optimiser. The IBM ILOG Cplex is freely available through Academic Initiative [here](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_|470|135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB). The [CBC](https://projects.coin-or.org/Cbc) solver is open source and freely available for any user. 

### References

[Terfve et al.](http://www.nature.com/articles/ncomms9033):

> Terfve, C. D. A., Wilkes, E. H., Casado, P., Cutillas, P. R., and Saez-Rodriguez, J. (2015). Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data. *Nature Communications*, 6:8033.

[Wilkes et al.](http://www.pnas.org/content/112/25/7719.abstract) (description of parts of the data)

> Wilkes, E. H., Terfve, C., Gribben, J. G., Saez-Rodriguez, J., and Cutillas, P. R. (2015). Empirical inference of circuitry and plasticity in a kinase signaling network. *Proceedings of the National Academy of Sciences of the United States of America,* 112(25):7719–24.

## Examples

### MTORi and PI3Ki-AKTi-MTORi Perturbation Analysis

Here we show the example from [Terfve et al.](http://www.nature.com/articles/ncomms9033) where PHONEMeS was used to model the perturbation effects of MTOR inhibition. 

We start first by loading the required packages and the necessary PHONEMeS inputs:

```R
# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)
library(dplyr)
library(readxl)
library(readr)
library(XML)

# Loading database and data-object
load(file = system.file("NetworKIN_noCSK_filt.RData", package="PHONEMeS"))
load(file = system.file("inputObj_Terfve.RData", package="PHONEMeS"))

```

Next we prepare the data objwcts which will be used as an input by PHONEMeS.

```R
# Preparation of the background network from the allD list of K/P-S interactions
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))

# Creating the list where we show the experimental conditions
conditions <- list(c("AKT1 - Control", "AKT2 - Control"), c("CAMK1 - Control", "CAMK2 - Control"),
                   c("EGFR1 - Control", "EGFR2 - Control"), c("ERK1 - Control", "ERK2 - Control"),
                   c("MEK1 - Control", "MEK2 - Control"), c("MTOR1 - Control", "MTOR2 - Control"),
                   c("P70S6K1 - Control", "P70S6K2 - Control"), c("PI3K1 - Control", "PI3K2 - Control"),
                   c("PKC1 - Control", "PKC2 - Control"), c("ROCK1 - Control", "ROCK2 - Control"))

names(conditions) <- c("AKT1_HUMAN", "KCC2D_HUMAN", "EGFR_HUMAN", "MK01_HUMAN",
                       "MP2K1_HUMAN", "MTOR_HUMAN",  "KS6B1_HUMAN", "PK3CA_HUMAN",
                       "KPCA_HUMAN", "ROCK1_HUMAN")

# For each experimental condition we assign a perturbation target
targets.P<-list(cond1=c("AKT1_HUMAN", "AKT2_HUMAN"), cond2=c("KCC2A_HUMAN", "KCC2B_HUMAN", "KCC2C_HUMAN", "KCC2D_HUMAN"), cond3=c("EGFR_HUMAN", "ERBB2_HUMAN"), 
                cond4=c("MK01_HUMAN", "MK03_HUMAN", "MK14_HUMAN"), cond5=c("MP2K1_HUMAN", "MP2K2_HUMAN"), cond6=c("MTOR_HUMAN"), cond7=c("KS6B1_HUMAN", "KS6B2_HUMAN"), 
                cond8=c("PK3CA_HUMAN", "PK3CD_HUMAN", "MTOR_HUMAN"), cond9=c("KPCA_HUMAN", "KPCB_HUMAN", "KPCG_HUMAN", "KPCE_HUMAN"), cond10=c("ROCK1_HUMAN", "ROCK2_HUMAN"))

```

Then finally we perform the PHONEMeS analysis for the MTOR inhibition experiment (condition 6).

**Attention:** Here the user must additionally provide the path to the executable CPLEX or CBC solver after having obtained the license and downloaded them. The path must be set on the ```solverPath``` variable (i.e. ```solverPath = "~/Documents/cplex"```, if the user wishes to use the *cplex* executable which he saved in the *Documents* direcotry).

```R
# Select experimental condition
experiments <- c(6) # for MTORi case

# Running PHONEMeS - cplex
# Run PHONEMeS with multiple solutions from CPLEX
resultsMulti <- runPHONEMeS(targets.P = targets.P, conditions = conditions, inputObj = inputObj, experiments = experiments, bg = bg, solver = "cplex", nSolutions = 100, nK = "no", solverPath = path_to_executable_solver)
write.table(x = resultsMulti, file = "MTORi_sif_cplex.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

nodesAttributes <- assignAttributes(sif = resultsMulti, dataGMM = inputObj, targets = targets.P[experiments], writeAttr = TRUE)

```

For running PHONEMeS by considering multiple conditions we can consider the case where we wish to retrieve consensus network solutions for when combining evidences from the PI3Ki-AKTi-MTORi inhibition experimental conditions. This analysis can be performed either by requesting multiple solutions directly from the CPLEX options or by performing PHONEMeS analysis multiple times by randomly downsampling the measurements and retrieving one solution for each iteration. These single solutions are then integrated into a combined network. With increasing number of experimental conditions, the number of possible solutions is expected to increase substantially and the Downsampling porcedure is advised to be used for these situations since it integrates single solutions in a more unbiased manner and the integrated network is sparser.

**Attention:** Here the user must additionally provide the path to the executable CPLEX or CBC solver after having obtained the license and downloaded them. The path must be set on the ```solverPath``` variable (i.e. ```solverPath = "~/Documents/cplex"```, if the user wishes to use the *cplex* executable which he saved in the *Documents* direcotry).

```R
# Select experimental condition
experiments <- c(1, 6, 8) # for PI3Ki-AKTi-MTORi case

# Running PHONEMeS - cplex
# Run PHONEMeS with multiple solutions from CPLEX
resultsMulti <- runPHONEMeS(targets.P = targets.P, conditions = conditions, inputObj = inputObj, experiments = experiments, bg = bg, solver = "cplex", nSolutions = 100, nK = "no", populate = 100, solverPath = path_to_executable_solver)
write.table(x = resultsMulti, file = "PI3Ki_AKTi_MTORi_sif_cplex.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
resultsMulti <- runPHONEMeS_Downsampling(targets.P = targets.P, conditions = conditions, inputObj = inputObj, experiments = experiments, bg = bg, nIter = 100, nK = "no")
write.table(x = resultsMulti, file = "PI3Ki_AKTi_MTORi_sif_downsampling.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

nodesAttributes <- assignAttributes(sif = resultsMulti, dataGMM = inputObj, targets = targets.P[experiments], writeAttr = TRUE)

```

### Time-Point Analysis for the Modelling of Endothelin Signalling

Here we show an example about how to perform PHONEMeS analsysis over phosphoproteomics time-course data. The example refers to the study from [Schaefer et al. 2019](https://www.embopress.org/doi/full/10.15252/msb.20198828). The PHONEMeS inputs used for this example have been prepared as described in detail in the separate [dedicated github repository](https://github.com/saezlab/EDN_phospho) of this study.

We start first by loading the required packages and the necessary PHONEMeS inputs for this example:

```R
# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)
library(dplyr)
library(readxl)
library(readr)
library(XML)

# Loading database and data-object
load(file = system.file("interactions_schaefer.RData", package="PHONEMeS"))
load(file = system.file("dataUACC257.RData", package="PHONEMeS"))

```

Next we prepare the data objwcts which will be used as an input by PHONEMeS.

```R
# Preparing background network as a PHONEMeS input
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataInput<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)

# Choose the conditions for each time-point
conditions <- list(c("tp_2min"), c("tp_10min"), c("tp_30min"), c("tp_60min"), c("tp_90min"))
names(conditions) <- c("tp_2min", "tp_10min", "tp_30min", "tp_60min", "tp_90min")

# Choose the targets for each time-point (EDNRB - the same)
targets.P <- list(tp_2min=c("EDNRB_HUMAN"), tp_10min=c("EDNRB_HUMAN"), tp_30min=c("EDNRB_HUMAN"), tp_60min=c("EDNRB_HUMAN"), tp_90min=c("EDNRB_HUMAN"))

# Next we assign the experimental conditions from which we get the measurements at each specific time-point
# In this case we have the same one experimental condition for each time-point
experiments <- list(tp1=c(1), tp2=c(2), tp3=c(3), tp4=c(4), tp5=c(5))

```

Next we perform the PHONEMeS analysis to obtain the time-course modelling of signalling. We perform this analysis 100 times where for each iteration we retain a random sample of measurements. Each solution at each iteration contains one time-course signalling model which in the end is then combined into one integrated network.

**Attention:** Here the user must additionally provide the path to the executable CPLEX or CBC solver after having obtained the license and downloaded them. The path must be set on the ```solverPath``` variable (i.e. ```solverPath = "~/Documents/cplex"```, if the user wishes to use the *cplex* executable which he saved in the *Documents* direcotry).

```R
# Running multiple time-point variant of PHONEMeS and retain only those interactions which have a weight higher than 20/appear at least 20 times in the
# separate solutions we have obtained out of the 100 runs we have set to perform (nIter=100).
set.seed(383789)
resultsMulti = runPHONEMeS_mult(targets.P = targets.P, conditions = conditions, inputObj = dataInput, experiments = experiments, bg = bg, nIter = 100)
nodeAttribudes <- assignAttributes(sif = resultsMulti[which(resultsMulti[, 2]>=20), ], dataGMM = dataInput, targets = targets.P, writeAttr = FALSE)

write.table(x = resultsMulti[which(resultsMulti[, 2]>=20), ], file = "ednrb_network.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = nodeAttribudes, file = "ednrb_attributes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
```

### PHONEMeS Upside-Down Analysis

As always, we start first by loading the required packages and the necessary PHONEMeS inputs for this example:
```R
# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)
library(dplyr)
library(readxl)
library(readr)
library(limma)
library(OmnipathR)
library(XML)

# Load data
load(file = system.file("colon_cancer_data.RData", package="PHONEMeS"))
load(file = system.file("colon_cancer_kinase_activities.RData", package="PHONEMeS"))

```

From the Colon Cancer data we perform a differential analysis to estimate differentially regulated phospho-sites between the tumour and health samples. For that we rely on the [Limma R-Package](https://bioconductor.org/packages/release/bioc/html/limma.html) and the built-in PHONEMeS functions specifically designed for this purpose.

```R
# Differential analysis with Limma

# Defining targets
targets <- matrix(data = , nrow = ncol(colon_cancer_data), ncol = 2)
targets[, 1] <- colnames(colon_cancer_data)
targets[, 2] <- sapply(strsplit(x = colnames(colon_cancer_data), split = "_", fixed = TRUE), "[[", 1)
colnames(targets) <- c("sample", "condition")
targets <- as.data.frame(targets)
targets$sample <- as.character(targets$sample)
targets$condition <- as.character(targets$condition)

# Performing the differential analysis
limmaRes <- runLimma(measurements = colon_cancer_data, targets = targets, comparisons = list(c(1, -2))) # tumor vs healthy
ttop <- topTable(fit = limmaRes[[1]], coef = 1, number = nrow(colon_cancer_data), adjust.method = "fdr")
ttop$ID <- rownames(ttop)

# Generating the table-top object from the Limma analysis result
ttopList <- list()
ttopList[[length(ttopList)+1]] <- ttop
names(ttopList) <- "tumour_vs_healthy"

```

Results from Limma can then be used to build the data-input object of PHONEMeS as below. Here we identify as significant those measurements with an
absolute *log2FC* threshold higher than 2 and adjusted p-value significance lower than 0.01.

```R
# Build PHONEMeS data-input object
dataInput <- buildInputs(tableTopList = ttopList, fcThresh = log2(2.5), pThresh = 0.05, idxID = 7, idxFC = 1, 
                         idxPval = 5, namesConditions = c("tumour_vs_healthy"))

```

PHONEMeS also needs a set Kinase/Phosphatase-Substrate to be used as a prior knowledge network-input which we train to the data. For that we can rely on OmniPath, and for PHONEMeS such a PKN can be built as follows (via the [OmnipathR R-package](https://github.com/saezlab/OmnipathR)):

```R
# Obtaining the PTMS relations online from OmnipathR package and retaining only phosphorylation/dephosphorylation modifications
ptms <- get_signed_ptms()
ptms <- ptms[which(ptms$modification%in%c("phosphorylation", "dephosphorylation")), ]

# Initializing an empty PKN for PHONEMEeS
bn <- matrix(data = , nrow = nrow(ptms), ncol = 8)
colnames(bn) <- c("S.AC", "S.ID", "K.AC", "K.ID", "res", "pos", "SID", "S.cc")

# Filling the initialized PKN with information from OmniPath
bn[, 1] <- ptms$substrate_genesymbol
bn[, 2] <- ptms$substrate_genesymbol
bn[, 3] <- ptms$enzyme_genesymbol
bn[, 4] <- ptms$enzyme_genesymbol
bn[, 5] <- ptms$residue_type
bn[, 6] <- ptms$residue_offset
bn[, 7] <- paste0("e", 1:nrow(bn))
bn[, 8] <- paste0(ptms$substrate_genesymbol, "_", ptms$residue_type, ptms$residue_offset)

# Making the PKN as a data-frame object
allD <- as.data.frame(bn)
allD$S.AC <- as.character(allD$S.AC)
allD$S.ID <- as.character(allD$S.ID)
allD$K.AC <- as.character(allD$K.AC)
allD$K.ID <- as.character(allD$K.ID)
allD$res <- as.character(allD$res)
allD$pos <- as.character(allD$pos)
allD$SID <- paste0("e", 1:nrow(allD))
allD$S.cc <- as.character(allD$S.cc)

```

Now preparing the other inputs for PHONEMeS and also selecting the most regulated kinases based on kinase enrichment analysis results (with absolute kinase activities higher than 12).

```R
# Identifying top most regulated kinases
targetKinases <- kinase_activities$kinase[which(abs(kinase_activities$activity)>=12)]

# Preparation of the background network from the allD list of K/P-S interactions
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))

# Creating the list where we show the experimental conditions
conditions <- list(c("tumour_vs_healthy"))

names(conditions) <- c("tumour_vs_healthy")

# For each experimental condition we assign a perturbation target
targets.P<-list(tumour_vs_healthy=targetKinases)

# Select experimental condition
experiments <- c(1)

```

Finally running the upside-down variant of PHONEMeS and saving the network results. The ```phonemes_ud``` output on this case will consist of a list of three objects: the downside network solution, the upside network solution as well as the combined network solution which integrates the two. These results can then be saved as matrices in a *txt* format which can then be used to visualize the network graph solutions.

**Attention:** Here the user must additionally provide the path to the executable CPLEX or CBC solver after having obtained the license and downloaded them. The path must be set on the ```solverPath``` variable (i.e. ```solverPath = "~/Documents/cplex"```, if the user wishes to use the *cplex* executable which he saved in the *Documents* direcotry).

```R
phonemes_ud <- runPHONEMeS_UD(targets.P = targets.P, conditions = conditions, dataGMM = dataInput, experiments = experiments, 
                              bg = bg, solverPath = path_to_executable_solver)

write.table(x = phonemes_ud$Downside, file = "phonemes_ud_down.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = phonemes_ud$Upside, file = "phonemes_ud_up.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = phonemes_ud$Combined, file = "phonemes_ud_combined.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

nodesAttributes <- assignAttributes(sif = phonemes_ud$Combined, dataGMM = dataInput, targets = targets.P)
write.table(x = nodesAttributes, file = "nodes_attributes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

```
