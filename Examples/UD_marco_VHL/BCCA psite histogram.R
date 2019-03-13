setwd("~/Documents/temp/ud_case_for_aurelien")

# Load packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(hash)
library(readr)
library(visNetwork)
library(threejs)
library(dplyr)
library(ggplot2)

study <- 4

uniprot_goa_metabolic_process <- as.data.frame(read_delim("~/Documents/PHONEMeS-ILP/Examples/VHL/uniprot_metabolic_process_GO.tab",
                                                          "\t", escape_double = FALSE, trim_ws = TRUE))
uniprot_goa_metabolic_process$`Gene names` <- gsub(" .*","",uniprot_goa_metabolic_process$`Gene names`)

uniprot_branched_chain_amino <- as.data.frame(read_delim("~/Documents/PHONEMeS-ILP/Examples/VHL/uniprot-_branched+chain+amino_-filtered-reviewed_yes+AND+organism__Hom--.tab",
                                                         "\t", escape_double = FALSE, trim_ws = TRUE))

uniprot_branched_chain_amino$`Gene names` <- gsub(" .*","",uniprot_branched_chain_amino$`Gene names`)

recon_store_genes_1 <- as.data.frame(read_delim("~/Documents/PHONEMeS-ILP/Examples/VHL/recon-store-genes-1.tsv",
                                                "\t", escape_double = FALSE, trim_ws = TRUE))

load("~/Documents/PHONEMeS-ILP/Examples/UD_marco_VHL/input_PHONEMES_UD.Rdata")

ttop_list <- list()
for( i in 1:length(for_enio))
{
  ttop <- for_enio[[i]][[2]]
  names(ttop)[1] <- "X1"
  ttop <- ttop[gsub("_.*","",ttop$X1) %in% uniprot_branched_chain_amino$`Gene names`,]
  ttop_list[[i]] <- ttop
}
ttop <- for_enio[[study]][[2]]
names(ttop)[1] <- "X1"
# ttop <- ttop[gsub("_.*","",ttop$X1) %in% recon_store_genes_1$symbol,]
# ttop <- ttop[gsub("_.*","",ttop$X1) %in% uniprot_goa_metabolic_process$`Gene names`,]
ttop <- ttop[gsub("_.*","",ttop$X1) %in% uniprot_branched_chain_amino$`Gene names`,]

t_vals <- as.data.frame(matrix(NA,0,3))
names(t_vals) <- c("X1","t","condition")

i <- 1
for (ttop in ttop_list)
{
  ttop$condition <- names(for_enio)[i]
  t_vals <- as.data.frame(rbind(t_vals,ttop[,c(1,4,8)]))
  i <- i+1
}

ggplot(t_vals, aes(x = X1, y = t, fill = condition)) + geom_bar(stat = "identity", position =  position_dodge()) + theme_minimal() + geom_hline(yintercept=range(2.2, -2.2), color='blue', size=1) + xlab("Phosphosite") + ylab("moderated t value") + ggtitle("PTMs deregulation in branched amino acid pathwway") +
  geom_text(aes( 3, 2.2, label = "FDR 0.1", vjust = -1), size = 5) + geom_text(aes( 3, -2.2, label = "FDR 0.1", vjust = 1.5), size = 5)
