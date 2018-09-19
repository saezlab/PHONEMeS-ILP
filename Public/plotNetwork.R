# Plot PHONEMeS networks using ggraph

library("dplyr")
library("igraph")
library("tidygraph")
library("ggraph")


plotNetwork <- function(resultsSIF, dataGMM, targets.P, ...){
  # list of nodes in the network to plot
  list_of_nodes <- unique(c(resultsSIF$Source, resultsSIF$Target))
  
  list_of_targets <- as.character(unlist(targets.P))
  
  # nodes that appear in the data
  GMM.ID <- dataGMM@IDmap
  sites <- intersect(GMM.ID$S.cc, list_of_nodes)

  # annotate nodes as perturbed (P) or targeted by the drug applied in the experiment (D)
  nodes_attributes <- data.frame(Species=list_of_nodes)
  nodes_attributes <- nodes_attributes %>% mutate(nodesP="")
  nodes_attributes <- nodes_attributes %>% mutate(nodesP=ifelse(Species %in% sites, "P", nodesP))
  nodes_attributes <- nodes_attributes %>% mutate(nodesP=ifelse(Species %in% list_of_targets, "D", nodesP))
  
  resultsSIF <- resultsSIF %>% mutate(from=Source, to=Target)
  tree <- tbl_graph(nodes=nodes_attributes, edges=resultsSIF, directed=TRUE)
  plt <- plotFormattedNetworkGgraph(tree)
  print(plt)
}

plotFormattedNetworkGgraph <- function(tree){
  # shape and colour definition linked to `nodesP`
  node_definition <- data.frame(nodesP=c("", "D", "P"), 
                                 shape=c("circle", "triangle", "diamond"),
                                 colour=c("grey", "red", "blue"))
  colour_def <- as.character(node_definition$colour)
  names(colour_def) <- as.character(node_definition$nodesP)
  shape_def <- as.character(node_definition$shape)
  names(shape_def) <- as.character(node_definition$nodesP)
  
  # plot network
  plt <- ggraph(tree, layout="tree") + 
    geom_node_point(size = 10, 
                    aes(shape=nodesP, colour=nodesP),
                    show.legend = FALSE) +
    scale_colour_manual(values = colour_def) +
    scale_shape_manual(values = shape_def) +
    geom_edge_link(arrow = arrow(length = unit(4, "mm")), 
                   start_cap = circle(3, "mm"),
                   end_cap = circle(3, "mm")) + 
    geom_node_text(aes(label=Species), repel=TRUE)
  
  return(plt)
}
