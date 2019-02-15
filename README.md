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
install.packages("PHONEMeS_0.2.7.tar.gz", repos=NULL)
```

The ILP formulation contained in this repository does not need installation, just download and run the desired scripts (see [Documentation](https://github.com/saezlab/PHONEMeS-ILP/tree/master/Examples) for more information).

Other supportive R packages needed are:

*igraph* which you can easily install by typing in R the below line:

```R
install.packages("igraph")
```

*XML* which can be downloaded [here](https://cran.r-project.org/src/contrib/XML_3.98-1.9.tar.gz) and then you can install by typing in R:

```R
install.packages("XML_3.98-1.9.tar.gz", repos=NULL)
```

### Usage

For a guide how to run a PHONEMeS analysis, please refer to the dedicated [documentation](https://github.com/saezlab/PHONEMeS-ILP/tree/master/Examples).

### CPLEX

PHONEMeS-ILP depends on the CPLEX solver, meaning that the user needs to obtain an IBM ILOG CPLEX license and then save the executable file on the working directory. The IBM ILOG CPLEX Optimization Studio license can be obtained for free by students and faculty staff/researchers [here](https://ibm.onthehub.com/) (under *Data & Analytics* -> *Software* -> *CPLEX*). Note that this downloads the whole CPLEX suite, but user only needs the executable file named `cplex`.

### References

[Terfve et al.](http://www.nature.com/articles/ncomms9033):

> Terfve, C. D. A., Wilkes, E. H., Casado, P., Cutillas, P. R., and Saez-Rodriguez, J. (2015). Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data. *Nature Communications*, 6:8033.

[Wilkes et al.](http://www.pnas.org/content/112/25/7719.abstract) (description of parts of the data)

> Wilkes, E. H., Terfve, C., Gribben, J. G., Saez-Rodriguez, J., and Cutillas, P. R. (2015). Empirical inference of circuitry and plasticity in a kinase signaling network. *Proceedings of the National Academy of Sciences of the United States of America,* 112(25):7719–24.
