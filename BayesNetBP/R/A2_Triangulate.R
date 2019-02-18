
Triangulate <- function(graph, elim.order){
  dag.graph <- igraph.from.graphNEL(graph, weight=FALSE)
  # iterate over all nodes by elimination order
  for (i in 1:length(elim.order)){
    node <- elim.order[i]
    neighbors <- names(neighbors(dag.graph, node, mode="all"))
    n.neighbors <- length(neighbors)
    if (n.neighbors>=2){
      # iterate over all pairs of neighbors
      for (k1 in 1:(n.neighbors-1)){
        for (k2 in (k1+1):n.neighbors){
          # if both of the neighbor appears later in EO & are not connected
          # then connect them
          if (which(elim.order==neighbors[k1])>i){
            if (which(elim.order==neighbors[k2])>i){
              if(!are_adjacent(dag.graph, neighbors[k1], neighbors[k2])){
                dag.graph <- add_edges(dag.graph, c(neighbors[k1],neighbors[k2]))
              }
            }
          }
        }
      }
    }
  }
  dag.tri <- igraph.to.graphNEL(dag.graph) 
  return(dag.tri)
}