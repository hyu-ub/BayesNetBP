#' @importFrom igraph get.edge.ids

Triangulate <- function(graph, elim.order){
  dag.graph <- igraph.from.graphNEL(graph, weight=FALSE)
  for(i in 1:length(elim.order)){
    node <- elim.order[i]
    neighbors <- names(neighbors(dag.graph, node, mode="all"))
    nn <- length(neighbors)
    if(nn >= 2){
      node_pairs <- generate_pairs(nn)
      for(j in 1:nrow(node_pairs)){
        n1 <- neighbors[node_pairs[j, 1]]
        n2 <- neighbors[node_pairs[j, 2]]
        if (which(elim.order==n1)>i && which(elim.order==n2)>i){
          if(get.edge.ids(dag.graph, c(n1, n2)) == 0){
            dag.graph <- add_edges(dag.graph, c(n1, n2))
          }
        }
      }
    }
  }
  dag.tri <- igraph.to.graphNEL(dag.graph)
  return(dag.tri)
}

# given a magnitude, generate all possible pairs of the numbers from 1 to magnitude inclusive. Order is not considered.
generate_pairs <- function(magnitude){
  if(magnitude <= 1) return(list())
  enumeration <- matrix(nrow = (magnitude-1)*magnitude/2, ncol = 2)
  sum = 0
  for(i in 1:(magnitude-1)){
    enumeration[i:(magnitude-1) + sum, 1] <- i
    enumeration[i:(magnitude-1) + sum, 2] <- (i+1):magnitude
    sum = sum + magnitude - i - 1
  }
  return(enumeration)
}
