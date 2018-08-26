Dataset from which the PHONEMeS objects were built (**dataGMM.RData** and **allD_MOUSE.RData** objects containing the scores assigne to each peptide across each condition and the background network to train).

The set of data in this directory contains:

**geneset_pdgf.txt** file containing the members of *BIOCARTA_PDGF_PATHWAY* from [MSigDB](http://software.broadinstitute.org/gsea/msigdb/cards/BIOCARTA_PDGF_PATHWAY.html)
  
**omnipath_mouse.txt** file containing the protein interactions from Omnipat for Mouse

**omnipath_mouse_enzyme-substrate_all.txt** file containing the k/p-s interactions for Mouse from Omnipath

**string_interactions_pdgfra.tsv** & **string_interactions_pdgfrb.tsv** files containing the direct proteins connected to PDGFR receptors from STRING

**ttop_list.RData** a table-top list containing the results from Limma analysis over the dataset provided in the paper.
