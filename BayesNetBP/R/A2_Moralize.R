#' @importFrom igraph as.undirected

Moralize <- function(graph){
  # graph <- emission1000$dag
  dag.graph <- igraph.from.graphNEL(graph, weight=FALSE)
  und.graph <- as.undirected(dag.graph, mode = "collapse")
  # iterate over all nodes 
  nds <- graph::nodes(graph)
  
  for (i in 1:length(nds)) {
    # i <- 1
    node <- nds[i]
    pa <- names(neighbors(dag.graph, node, mode="in"))
    n.pa <- length(pa)
    
    if (n.pa>=2) {
      for (k1 in 1:(n.pa - 1)) {
        for (k2 in (k1+1):n.pa) {
          if(!are_adjacent(und.graph, pa[k1], pa[k2])){
            und.graph <- add_edges(und.graph, c(pa[k1],pa[k2]))
          }
        }
      }
    }
  }
  
  nel.mor <- igraph.to.graphNEL(und.graph) 
  # x11(); plot(nel.mor)
  return(nel.mor)
  
}
