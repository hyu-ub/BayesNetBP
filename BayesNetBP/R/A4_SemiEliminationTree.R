#' @importFrom igraph igraph.from.graphNEL
SemiEliminationTree <- function(graph, cluster.sets, node.class, elim.order){
  cs.graph <- igraph.from.graphNEL(graph)
  i <- 1

  while (i < length(cluster.sets)){
    # cat("one iteration")
    cs.names <- names(cluster.sets)
    cluster <- cluster.sets[[i]]

    this.name <- cs.names[i]
    this.parent <- neighbors(cs.graph, v=this.name, mode="in")$name

    find_cluster <- FALSE
    is.discrete <- prod(node.class[cluster]) # if this cluster is dicrete (all members are discete)
    if (is.discrete){
      # if it's discrete, search for the following clusters
      for(j in (i+1):length(cluster.sets)){
        if (prod(node.class[cluster.sets[[j]]])){
          # if jth cluster is discrete
          if (all(cluster %in% cluster.sets[[j]])){
            # and it contains this cluster, set the indicator true
            find_cluster <-  TRUE
            break
            # stop searching
          }
        }
      }
    }

    # if can not found such a cluster j, then proceed to next iteration
    if(!find_cluster){i <- i+1; next}

    # if such cluster j is found, do the following:

    # find all the children of this cluster in the elimination tree
    cluster.children <- neighbors(cs.graph, v=cs.names[i], mode="out")$name

    # find a such a child among them,
    # whose elimination node should appear as late as possible in the order
    selected.child <- 0
    selected.child.name <- NULL
    current <- 0

    for(j in 1:length(cluster.children)){
      this.child <- cluster.children[[j]]
      child.mem <- cluster.sets[[this.child]]
      is.discrete <- prod(node.class[child.mem])
      if (is.discrete & all(cluster %in% child.mem)) {
        this.scan <- which(elim.order == cs.names[i])

        # if the elimination node of this child appears later
        # then update the selected child
        if(this.scan > current){
          current <- this.scan ##
          selected.child <- j;
          selected.child.name <- cluster.children[j]
        }
      }
    }

    # merge this cluster into its selected child by
    # 1. deleting this cluster
    cs.graph <- delete_vertices(cs.graph, this.name)

    # 2. update parents and children of involved nodes
    if(length(this.parent)!=0){ # if this cluster has a parent
      for (j in 1:length(this.child)){
        if(this.child!=selected.child.name){
          # let the parent of all its child be its selected child, except for the selected child itself
          cs.graph <- add_edges(cs.graph, c(selected.child.name, this.child[j]))
        } else {
          # let its parent be the parent of its selected child
          cs.graph <- add_edges(cs.graph, c(this.parent, this.child[j]))
        }
      }
    } else { # if it does not have a parent, then update as usual,
      # but no need to update the parent of its selected child
      for (j in 1:length(cluster.children)){
        if(cluster.children[j]!=selected.child.name){
          cs.graph <- add_edges(cs.graph, c(selected.child.name, cluster.children[j]))
        }
      }
    }

    # update the name vector of the tree
    selected.child.index <- which(cs.names==selected.child.name)
    cluster.sets[this.name] <- cluster.sets[selected.child.index]
    cs.names[which(cs.names==this.name)] <- selected.child.name

    cluster.sets <- cluster.sets[-selected.child.index]
    cs.names <- cs.names[-selected.child.index]
    names(cluster.sets) <- cs.names

  }

  # V(cs.graph)$name <- unlist(lapply(cluster.sets, paste0, collapse=","))
  # x11(); plot(igraph.to.graphNEL(cs.graph))

  return(igraph.to.graphNEL(cs.graph))

}
