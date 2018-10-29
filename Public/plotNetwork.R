# Plot PHONEMeS networks using ggraph

library("dplyr")
library("igraph")
library("tidygraph")
library("ggraph")


plotNetwork <- function(resultsSIF, dataGMM, targets.P, print=TRUE, ...){
  # annotate nodes as perturbed (P) or targeted by the drug applied in the experiment (D)
  nodes_attributes <- annotate_nodes_from_resultsSIF(resultsSIF, dataGMM, targets.P)
    
  # use `tidygraph` for representing the network
  resultsSIF <- resultsSIF %>% mutate(from=Source, to=Target)
  tree <- tbl_graph(nodes=nodes_attributes, edges=resultsSIF, directed=TRUE)
  plt <- plotFormattedNetworkGgraph(tree, ...)
  
  # plot directly or return the ggraph plot object, which can be modified further
  if (print){
    print(plt)
  } else {
    return(plt)
  }
}

annotate_nodes_from_resultsSIF <- function(resultsSIF, dataGMM, targets.P){
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
  
  return(nodes_attributes)
}

plotFormattedNetworkGgraph <- function(tree, layout="tree", repel=TRUE){
  # shape and colour definition linked to `nodesP`
  node_definition <- data.frame(nodesP=c("", "D", "P"), 
                                 shape=c("circle", "triangle", "diamond"),
                                 colour=c("grey", "red", "blue"))
  colour_def <- as.character(node_definition$colour)
  names(colour_def) <- as.character(node_definition$nodesP)
  shape_def <- as.character(node_definition$shape)
  names(shape_def) <- as.character(node_definition$nodesP)
  
  # plot network
  plt <- ggraph(tree, layout=layout) + 
    geom_node_point(size = 10, 
                    aes(shape=nodesP, colour=nodesP),
                    show.legend = FALSE) +
    scale_colour_manual(values = colour_def) +
    scale_shape_manual(values = shape_def) +
    geom_edge_link(arrow = arrow(length = unit(4, "mm")), 
                   start_cap = circle(3, "mm"),
                   end_cap = circle(3, "mm")) + 
    #geom_node_text(aes(label=Species), repel=repel) +
    geom_node_label(aes(label=Species), repel=repel) +
    theme_graph()
  
  return(plt)
}

# # test code for the development of the ggraph plotting function
# tree <- create_tree(20, 3)
# tree <- tree %>%
#   activate(nodes) %>%
#   mutate(Species=sort(unique(c(.E()$from, .E()$to)))) %>%
#   mutate(nodesP="") %>%
#   mutate(nodesP=ifelse(Species>15, "D", nodesP)) %>%
#   mutate(nodesP=ifelse(Species==1, "P", nodesP))
# tree %>% plotFormattedNetworkGgraph() %>% print()
# # if the tree contains a feedback loop, a warning message is displayed
# # and the final plot may not be usable
# tree2 <- tree %>% bind_edges(data.frame(from=20, to=2))

# TODO: Finish code for plotting with Cytoscape
# note: the use of yFiles layouts is not possible because of licensing issues
# See https://github.com/cytoscape/cyREST/issues/36
# Do not execute this part of the code as long as the implementation using RCy3
# is not finished.
if (FALSE){
  # Plotting PHONEMeS networks using RCy3
  library("RCy3")
  
  plotNetworkRCy3 <- function(resultsSIF, dataGMM, targets.P, ...){
    
    is_cytoscape_running()
    
    # annotate nodes as perturbed (P) or targeted by the drug applied in the experiment (D)
    nodes_attributes <- annotate_nodes_from_resultsSIF(resultsSIF, dataGMM, targets.P)
    
    # use `tidygraph` for representing the network
    resultsSIF <- resultsSIF %>% mutate(from=Source, to=Target)
    tree <- tbl_graph(nodes=nodes_attributes, edges=resultsSIF, directed=TRUE)
    session_id <- plotFormattedNetworkRCy3(tree, ...)
    
    
    # # list of nodes in the network to plot
    # list_of_nodes <- unique(c(resultsSIF$Source, resultsSIF$Target))
    # 
    # list_of_targets <- as.character(unlist(targets.P))
    # 
    # # nodes that appear in the data
    # GMM.ID <- dataGMM@IDmap
    # sites <- intersect(GMM.ID$S.cc, list_of_nodes)
    # 
    # # annotate nodes as perturbed (P) or targeted by the drug applied in the experiment (D)
    # nodes_attributes <- data.frame(Species=list_of_nodes)
    # nodes_attributes <- nodes_attributes %>% mutate(nodesP="")
    # nodes_attributes <- nodes_attributes %>% mutate(nodesP=ifelse(Species %in% sites, "P", nodesP))
    # nodes_attributes <- nodes_attributes %>% mutate(nodesP=ifelse(Species %in% list_of_targets, "D", nodesP))
    # 
    # plotFormattedNetwork(network_df, ...)
  }
  
  
  plotFormattedNetworkRCy3 <- function(network_df, title="network", collection="Example"){
    network_df <- network_df %>% activate(nodes) %>% mutate(id=Species)
    network_df <- network_df %>% activate(edges) %>% 
      mutate(source=Source, target=Target, interaction=Interaction)
    #session_id <- createNetworkFromIgraph(network_df, title=title, collection=collection)
    node_df <- network_df %>% activate(nodes) %>% as.data.frame()
    edge_df <- network_df %>% activate(edges) %>% as.data.frame()
    session_id <- createNetworkFromDataFrames(node_df, edge_df, title=title, collection=collection)
    return(session_id)
  }
  
  is_cytoscape_running <- function(){
    tryCatch(
      {
        cytoscapePing()
        return(TRUE)
      }, error=function(error_message){
        message("Cytoscape needs to be running in the background.")
        message("Please, open Cytoscape before calling this function.")
        message(error_message)
        return(FALSE)
      }
    )
  }
  
  applyPhonemesStyle() <- function(){
    
    # node shape
    column <- 'nodesP'
    values <- c ('',  'D','P')
    shapes <- c ('ELLIPSE', 'DIAMOND', 'HEXAGON')
    setNodeShapeMapping (column, values, shapes)
    
    # node color
  }
}
