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
# for( i in 1:length(for_enio))
# {
#   ttop <- for_enio[[i]][[2]]
#   names(ttop)[1] <- "X1"
#   ttop_list[[i]] <- ttop
# }
ttop <- for_enio[[study]][[2]]
names(ttop)[1] <- "X1"
# ttop <- ttop[gsub("_.*","",ttop$X1) %in% recon_store_genes_1$symbol,]
# ttop <- ttop[gsub("_.*","",ttop$X1) %in% uniprot_goa_metabolic_process$`Gene names`,]
ttop <- ttop[gsub("_.*","",ttop$X1) %in% uniprot_branched_chain_amino$`Gene names`,]

ttop_list[[1]] <- ttop
names(ttop_list) <- c("cond1")

mappingTable <- as.data.frame(matrix(NA,length(ttop_list[[1]][,1]),6))
names(mappingTable) <- c("dataID","UPID","site","res","pos","S.cc")

mappingTable$dataID <- ttop_list[[1]][,1]
mappingTable$UPID <- gsub("_.*","",mappingTable$dataID)
mappingTable$site <- gsub(".*_","",mappingTable$dataID)
mappingTable$res <- gsub("[0-9]*","",gsub(".*_","",mappingTable$dataID))
mappingTable$pos <- mappingTable$pos <- gsub("[A-Za-z]","",mappingTable$site)
mappingTable$S.cc <- mappingTable$dataID

source("buildInputs.R")
ttop_list[[1]] <- as_tibble(ttop_list[[1]])
dataInput <- buildInputs(tableTopList = ttop_list, organism = "HUMAN", pThresh = 0.05, mappingTable = mappingTable)

pkn <- for_enio[[study]][[1]]

allD <- as.data.frame(matrix(NA,length(for_enio[[study]][[1]][,1]),8))
names(allD) <- c("S.AC","S.ID","K.AC","K.ID","res","pos","SID","S.cc")

allD$S.AC <- gsub("_.*","",pkn$substrate_genesymbol)
allD$S.ID <- allD$S.AC
allD$K.AC <- as.character(pkn$enzyme_genesymbol)
allD$K.ID <- allD$K.AC
allD$res <- gsub("[0-9].*","",gsub(".*_","",pkn$substrate_genesymbol))
allD$pos <- gsub("[A-Za-z]*","",gsub(".*_","",pkn$substrate_genesymbol))
allD$SID <- paste("e",(1:length(allD[,1])),sep = "")
allD$S.cc <- pkn$substrate_genesymbol

allD <- allD[complete.cases(allD), ]
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=dataInput$res, IDmap=dataInput$IDmap, resFC=dataInput$resFC)

kinases <- for_enio[[study]][[3]]
# kinases <- kinases[kinases %in% uniprot_goa_metabolic_process$`Gene names`]
targets.P <- list(cond1=kinases)

conditions <- list(c("cond1"))

names(conditions) <- "cond1"

data.P <- dataBycond(dataGMM, bg, scaled = TRUE, rowBycond = conditions)
experiments <- conditions

show(data.P)

speciesP(data.P)




source("buildDataMatrix.R")
source("ilpFunctions.R")
source("PKN_list.R")
source("build_Nw_Inv.R")
source("build_PKN_Inv.R")
source("runPHONEMeS.R")
source("runPHONEMeS-UD.R")
source("build_Nw.R")
source("build_PKN.R")

resUD <- runPHONEMeSUD(dataInput = dataInput, bg = bg, targets.P = targets.P, data.P = data.P, experiments = experiments)
# resUD <- runPHONEMeS(dataInput = dataInput, bg = bg, targets.P = targets.P, data.P = data.P, experiments = experiments)


###############


comb_metastasis_vs_cancer <- as.data.frame(resUD$Combined)

names(comb_metastasis_vs_cancer) <- c("from","sign","to")
comb_metastasis_vs_cancer$arrows <- "to"

attributes_metastasis_vs_cancer <- as.data.frame(resUD$nodesAttributes)
names(attributes_metastasis_vs_cancer) <- c("id","state")

attributes_metastasis_vs_cancer$shape <- "circle"
attributes_metastasis_vs_cancer[grepl("_",attributes_metastasis_vs_cancer$id),3] <- "triangle"

attributes_metastasis_vs_cancer$color <- "grey"
attributes_metastasis_vs_cancer[attributes_metastasis_vs_cancer$shape == "circle" & !is.na(attributes_metastasis_vs_cancer$state),4] <- "green"
attributes_metastasis_vs_cancer[attributes_metastasis_vs_cancer$shape == "triangle" & !is.na(attributes_metastasis_vs_cancer$state),4] <- "blue"

attributes_metastasis_vs_cancer$size <- 10
attributes_metastasis_vs_cancer[attributes_metastasis_vs_cancer$shape == "triangle",5] <- 10

attributes_metastasis_vs_cancer$label <- attributes_metastasis_vs_cancer$id

attributes_metastasis_vs_cancer <- attributes_metastasis_vs_cancer[,-2]

visNetwork(nodes = attributes_metastasis_vs_cancer, edges = comb_metastasis_vs_cancer) 
# visNetwork(nodes = attributes_metastasis_vs_cancer, edges = comb_metastasis_vs_cancer) %>%
#   visIgraphLayout(layout = "layout_as_tree")