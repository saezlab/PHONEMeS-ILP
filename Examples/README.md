# Running PHONEMeS

Running PHONEMeS-ILP is easy and straightforward once you have downloaded and installed [R](https://www.r-project.org/) on your computer and get a license for the [IBM CPLEX Optimizer](https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) which can be obtained for free for academic use. Also required is to have installed the supportive [R packages](https://github.com/saezlab/PHONEMeS-ILP) mentioned.

The pipeline is the same for all the case studies present in the Examples repository:

## Data input

* Download and save the data in a local folder where you wish to perform the analysis. The data folder will typicaly contain a background network containing a set of possible kinase to substrate interactions which we might obtain from multiple online data-bases such as Omnipath. Also download and save the files present on the Script folder.

* By running the ```buildDataMatrix.R```, ```ilpFunctions.R``` and ```ilpFunctions2.R``` scripts you are calling the necessay functions needed to write the objective function, bounds and constraints of the ILP problem.

* Then by running the ```executionScript.R``` you execute all the functions called from which you obtain the final *resultsSIF.txt* file which is basically the result that PHONEMeS gives you and which represents the optimal reconstructed pathway activity model.
