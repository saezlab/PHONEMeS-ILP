# PHONEMeS-ILP
ILP implementation of PHONEMeS - Enio GJERGA

## PHONEMeS Overview
PHONEMeS (PHOsphorylation NEtworks for Mass Spectrometry) is a method to model signaling networks based on untargeted phosphoproteomics mass spectrometry data and kinase/phosphatase-substrate interactions.

This package contains the R package and accompanying scripts that implement the method as well as several examples of how to run a PHONEMeS analysis.

The input for PHONEMeS consists of phosphoproteomic data after treatment with kinase inhibitors. Several approaches, such as Gaussian Mixture Modeling, can then be used to find phosphosites that exhibit a naturally Boolean behaviour with two populations, representing a control and a perturbed state. The data are mapped unto a kinase/phosphatase-substrate network taken from several dedicated databases such as Omnipath. PHONEMeS-ILP then optimizes the network and extracts possible paths connecting inhibited kinases and perturbed phosphosites through this new ILP (Integer Linear Programming) implementation.
