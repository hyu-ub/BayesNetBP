#' Convert a bn object to graphNEL object
#'
#' Convert a bn object to graphNEL object while removing isolated nodes
#'
#' @param graph_bn a \code{bn} object of Bayesian network
#' @return a \code{graphNEL} object
#' 
#' @author Han Yu
#'
#' @importFrom igraph igraph.from.graphNEL igraph.to.graphNEL induced_subgraph degree
#' @importFrom bnlearn as.graphNEL
#' 
#' @export

bn_to_graphNEL <- function(graph_bn) {
  
  graph.graphNEL <- bnlearn::as.graphNEL(graph_bn)
  graph.igraph <- igraph.from.graphNEL(graph.graphNEL)
  
  # remove isolated nodes
  deg <- igraph::degree(graph.igraph)
  nodes <- names(deg)[deg>0] 
  graph.igraph.sub <- induced_subgraph(graph.igraph, nodes)
  graph.graphNEL.sub <- igraph.to.graphNEL(graph.igraph.sub)
  
  return(graph.graphNEL.sub)
  
}



