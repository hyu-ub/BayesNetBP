#' @importFrom igraph as.undirected
#' @importFrom graph inEdges
#' @importFrom igraph simplify
Moralize <- function(graph){

  dag_nodes <- nodes(graph)

  und.graph <- as.undirected(igraph.from.graphNEL(graph, weight=FALSE), mode = "collapse")

  for(i in 1:length(dag_nodes)){
    parents <- inEdges(dag_nodes[i], graph)[[dag_nodes[i]]]
    parents_length <- length(parents)
    if(parents_length >= 2){
      for(p1 in 1:(parents_length-1)){
        for(p2 in (p1+1):(parents_length)){

            und.graph <- add_edges(und.graph, c(parents[p1], parents[p2])) # are_adjacent is expensive, so we don't use that.

        }
      }
    }

  }
  und.graph <- simplify(und.graph, remove.loops=FALSE)
  nel.mor <- igraph.to.graphNEL(und.graph)
  return(nel.mor)

}
