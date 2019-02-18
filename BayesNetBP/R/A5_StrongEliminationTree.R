
StrongEliminationTree <- function(cs, elim.order){
  
  # a vector to store the parent of each cluster
  parentnodes <- c()

  
  # iterate over all clusters
  for (i in length(cs):1){
    cluster <- cs[[i]]
    
    elim.ind <- which(elim.order == names(cs)[i])
    
    # Find a member such that
    # (1) its order is as small as possible, 
    # (2) but it appears after the elim node of this cluster in order
    # after finding this member, find the cluster (B) with this member as its elimination node
    # set cluster B as the parent of cluster A
    
    eo.pos <- match(cluster, elim.order)
    names(eo.pos) <- cluster
    eo.pos <- eo.pos[eo.pos>elim.ind]
    if (length(eo.pos)>0) {
      parentnodes[i] <- names(eo.pos)[which.min(eo.pos)]
    }
    
    
  }
  
  # construct a graphNEL object of the strong elimination tree
  
  names(parentnodes) <- names(cs)
  
  nodes <- names(parentnodes)
  Adj <- matrix(0, length(nodes), length(nodes))
  colnames(Adj) <- nodes
  rownames(Adj) <- nodes
  
  cs.graph <- graph_from_adjacency_matrix(Adj, mode = "directed")
  for (i in 1:length(nodes)){
    if(!is.na(parentnodes[i])){
      cs.graph <- add_edges(cs.graph, c(parentnodes[i], nodes[i]))
    }
  }
  
  cs.graph <- igraph.to.graphNEL(cs.graph)
  
  return(cs.graph)
}