#' @importFrom igraph as.undirected
#' @importFrom graph inEdges
#' @importFrom igraph simplify

Moralize <- function(graph){
  dag_nodes <- nodes(graph)
  und.graph <- as.undirected(igraph.from.graphNEL(graph, weight=FALSE), mode = "collapse")
  for(i in 1:length(dag_nodes)){
    parents <- inEdges(dag_nodes[i], graph)[[dag_nodes[i]]]
    np <- length(parents)
    if(np<2) next
    for(j in 1:(np-1)){
      for(k in (j+1):(np)){
          und.graph <- add_edges(und.graph, c(parents[j], parents[k])) 
      }
    }
  }
  und.graph <- simplify(und.graph, remove.loops=FALSE)
  nel.mor <- igraph.to.graphNEL(und.graph)
  return(nel.mor)
}
