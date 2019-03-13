library(readr)
library(visNetwork)
library(threejs)
library(omicToolsTest)
library(piano)
library(ggplot2)

string_interactions <- as.data.frame(read_delim("Documents/PHONEMeS-ILP/Examples/UD_marco_VHL/string_interactions.tsv", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE))

nodes <- data.frame(unique(c(string_interactions$`#node1`, string_interactions$node2)))
names(nodes) <- "id"
nodes$label <- nodes$id
nodes$level <- 3
nodes[nodes$id == "PPARGC1A","level"] <- 1
nodes[nodes$id %in% c("PPARG","ESRRA"),"level"] <- 2


edges <- string_interactions[,c(1,2)]
names(edges) <- c("from","to")

visNetwork(nodes = nodes, edges = edges) %>%
  visHierarchicalLayout()

pathways <- gmt_to_csv("~/Documents/PHONEMeS-ILP/Examples/UD_marco_VHL/c2.cp.v6.2.symbols.gmt")
pathways <- pathways[!grepl("PID",pathways$term),]

pathway_GSC <- loadGSC(pathways)

pathway_hyper <- runGSAhyper(genes = nodes$id, universe = as.character(unique(pathways$gene)), gsc = pathway_GSC)
pathway_hyper_df <- as.data.frame(pathway_hyper$pvalues)
names(pathway_hyper_df) <- "-log10(pvalue)"
pathway_hyper_df$`-log10(pvalue)` <- -log10(pathway_hyper_df$`-log10(pvalue)`)
pathway_hyper_df$pathway <- gsub("[_]"," ",row.names(pathway_hyper_df))     
pathway_hyper_df <- pathway_hyper_df[order(pathway_hyper_df$`-log10(pvalue)`, decreasing = T),]
pathway_hyper_df$pathway <- factor(pathway_hyper_df$pathway, levels = pathway_hyper_df$pathway)
pathway_hyper_df_top <- pathway_hyper_df[c(1:15),]

ggplot(pathway_hyper_df_top, aes(x = pathway, y = `-log10(pvalue)`)) + geom_bar(stat = "identity", fill = "royalblue4") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.margin = unit(c(1,1,1,10), "cm"))

write_csv(pathway_hyper_df[,c(2,1)], "~/Documents/PHONEMeS-ILP/Examples/UD_marco_VHL/pathway_hyper.csv")
