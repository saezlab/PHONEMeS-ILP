# PHONEMeS-ILP
ILP implementation of PHONEMeS - Enio GJERGA

**PHONEMeS** (**PHO**sphorylation **NE**tworks for **M**ass **S**pectrometry) is a method to model signalling networks based on untargeted phosphoproteomics mass spectrometry data and kinase/phosphatase-substrate interactions. 
Please see [Terfve et al.](http://www.nature.com/articles/ncomms9033) for an explanation of the methodolgy.

This repository contains the scripts of the ILP (Integer Linear Programming) implementation of the [PHONEMeS R package](https://github.com/saezlab/PHONEMeS/tree/master/Package) and accompanying scripts that implement the method. ILP is a mathematical optimisation problem in which the objective function and constraints are linear, while the variables are integers.

### License

Distributed under the GNU GPLv2 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/PHONEMeS/blob/master/LICENSE.txt) or copy at https://www.gnu.org/licenses/gpl-2.0.html.

### Installation

Before using the method, please install the current R package for PHONEMeS. For installation, download the tar file of the package and type in R:

```R
install.packages("PHONEMeS_0.2.7.tar.gz", repos=NULL)
```

Other supportive R packages needed are:

*igraph* which you can easily install by typing in R the below line:
```R
install.packages("igraph")
```

*BioNet* which you can easily install by typing in R the below line:
```R
source("https://bioconductor.org/biocLite.R")
biocLite("BioNet")
```

*XML* which can be downloaded [here](https://cran.r-project.org/src/contrib/XML_3.98-1.9.tar.gz) and then you can install by typing in R:

```R
install.packages("XML_3.98-1.9.tar.gz", repos=NULL)
```

### Usage

For a guide how to run a PHONEMeS analysis, please refer to the dedicated [documentation](https://github.com/saezlab/PHONEMeS-ILP/tree/master/Examples).

### CPLEX

PHONEMeS-ILP is CPLEX depndent meaning that the user needs to obtain an IBM ILOG CPLEX licence and then save the executable file to the working directory. The IBM ILOG CPLEX Optimization Studio license can be obtained for free by:

* Students: https://ibm.onthehub.com/WebStore/OfferingDetails.aspx?o=9b4eadea-9776-e611-9421-b8ca3a5db7a1

* Teachers, researchers and university staff: https://ibm.onthehub.com/WebStore/OfferingDetails.aspx?o=6fcc1096-7169-e611-9420-b8ca3a5db7a1

### References

[Terfve et al.](http://www.nature.com/articles/ncomms9033):

> Terfve, C. D. A., Wilkes, E. H., Casado, P., Cutillas, P. R., and Saez-Rodriguez, J. (2015). Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data. *Nature Communications*, 6:8033.

[Wilkes et al.](http://www.pnas.org/content/112/25/7719.abstract) (description of parts of the data)

> Wilkes, E. H., Terfve, C., Gribben, J. G., Saez-Rodriguez, J., and Cutillas, P. R. (2015). Empirical inference of circuitry and plasticity in a kinase signaling network. *Proceedings of the National Academy of Sciences of the United States of America,* 112(25):7719–24.
