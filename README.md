# PHONEMeS-ILP
ILP implementation of PHONEMeS - Enio GJERGA

## PHONEMeS Overview
PHONEMeS (PHOsphorylation NEtworks for Mass Spectrometry) is a method to model signaling networks based on untargeted phosphoproteomics mass spectrometry data and kinase/phosphatase-substrate interactions.

This package contains the R package and accompanying scripts that implement the method as well as several examples of how to run a PHONEMeS analysis.

The input for PHONEMeS consists of phosphoproteomic data after treatment with kinase inhibitors. Several approaches, such as Gaussian Mixture Modeling, can then be used to find phosphosites that exhibit a naturally Boolean behaviour with two populations, representing a control and a perturbed state. The data are mapped unto a kinase/phosphatase-substrate network taken from several dedicated databases such as Omnipath. PHONEMeS-ILP then optimizes the network and extracts possible paths connecting inhibited kinases and perturbed phosphosites through this new ILP (Integer Linear Programming) implementation.

Form more information about PHONEMeS, pleas visit it's official dedicated web-page: https://saezlab.github.io/PHONEMeS/

Main References:
Terfve, C. D. A., Wilkes, E. H., Casado, P., Cutillas, P. R., and Saez-Rodriguez, J. (2015). [Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data](https://www.nature.com/articles/ncomms9033). Nature Communications, 6:8033.

Wilkes, E. H., Terfve, C., Gribben, J. G., Saez-Rodriguez, J., and Cutillas, P. R. (2015). [Empirical inference of circuitry and plasticity in a kinase signaling network](http://www.pnas.org/content/112/25/7719.abstract). Proceedings of the National Academy of Sciences of the United States of America, 112(25):7719–24.


## ILP Implementation
