###########################################
## Propagate
###########################################

# targets is a vector of ANY nodes in the graph
propagate.worker <- function(tree.graph, potentials, cluster.sets, targets = NA){

  decomposed_tree <- igraph::decompose(tree.graph)
  if(any(!is.na(targets))){ # if any elements in targets are not NA
    cleaned_targets <- targets[!is.na(targets)] # remove NA
    cleaned_targets <- cleaned_targets[cleaned_targets %in% unlist(cluster.sets)] # remove invalid targets

    temp <- length(decomposed_tree)
    for(j in 1:temp){ # for every tree component
      if(!any( cleaned_targets %in% names(V(decomposed_tree[[j]])) )){ # if NONE of the cleaned targets are in the component
        decomposed_tree[[j]] <- NA # set it to NA for later removal
      }
    }

    decomposed_tree <- decomposed_tree[!is.na(decomposed_tree)] # remove all NA components(from above)
  }

  component_number <- length(decomposed_tree)
  worker_results <- vector("list", component_number)
  if(component_number > 0){
    for(i in 1:component_number){
      temp_names <- names(V(decomposed_tree[[i]]))

      # we know propagate.worker_orig works for single-component graphs
      # can exploit this and propagate.worker_orig on each component and then combine the results
      worker_results[[i]] <- propagate.worker_orig(decomposed_tree[[i]], potentials[temp_names], cluster.sets[temp_names])
    }
    # print(worker_results)
    return(unlist(worker_results, recursive = FALSE))
  }
  else{
    return(propagate.worker(tree.graph, potentials, cluster.sets))
  }
}

propagate.worker_orig <- function(tree.graph, potentials, cluster.sets){

  # tree.graph <- tree.sub.graph; potentials <- potentials.sub; cluster.sets <- discrete.sets

  cluster.tree <- list(
    # bn=dag,
    tree=tree.graph,
    clusters=cluster.sets,
    # assignment=asgn,
    collected=c(), active=c(), potentials=potentials, joint=potentials)

  clusters <- names(potentials)
  result <- list()

  ## NEW version of getting joints

   # print(clusters[1])

  # collect
  ce <- CollectEvidence(cluster.tree, clusters[1])
  # reset active nodes
  ce$active <- c()
  # distribute
  de <- DistributeEvidence(ce, clusters[1])
  result <- de$joint

  return(result)
}

# ppgt <- Propagate(tree.sub.graph, potentials.sub, discrete.sets)

###################################################

Absorb <- function(absorbedTo, absorbedFrom, separator, distribute=FALSE){
  pot1 <- absorbedFrom
  pot2 <- absorbedTo

  # pot2 <- cluster.tree$potentials[["Cyp2b10"]]; pot1 <- cluster.tree$potentials[["HDL"]];
  # separator <- intersect( cluster.tree$clusters[["Cyp2b10"]],  cluster.tree$clusters[["HDL"]])

  inter.var <- separator
  sep <- marginalize.discrete(pot1, inter.var)

  results <- list()
  if (distribute) {
    results[[1]] <- NULL
  } else {
    results[[1]] <- factor.divide(pot1, sep)
  }

  results[[2]] <- factor.product(pot2, sep)
  return(results)
}

###########################################
## Collect evidence
###########################################

CollectEvidence <- function(cluster.tree, node){
  process_queue <- c(node) # nodes to be processed - end nodes are last (back)
  node_stack <- c(list(c(node, NA))) # nodes to be absorbed - end nodes are first (front)
  # add in pairs: child, parent pairs

  processed <- c(node)

  while(length(process_queue) > 0){
    neighbors <- names(neighbors(cluster.tree$tree, process_queue[1], mode = "all")) # get neighbors
    unprocessed_neighbors <- neighbors[!(neighbors %in% processed)] # cull for neighbors not in stack

    if(length(unprocessed_neighbors) > 0){
      temp <- cbind(unprocessed_neighbors, rep(process_queue[1], length(unprocessed_neighbors)))
      node_stack <- c(split(temp, 1:nrow(temp)), node_stack) # add neighbors to stack top

      process_queue <- c(process_queue[-1], unprocessed_neighbors) # remove first node and add neighbors to queue back
    }
    else{
      process_queue <- process_queue[-1]
    }

    processed <- c(processed, unprocessed_neighbors)

  }

  for(i in 1:(length(node_stack)-1)){
    child <- node_stack[[i]][1]
    parent <- node_stack[[i]][2]


    abb <- Absorb(cluster.tree$potentials[[parent]], cluster.tree$potentials[[child]],
                  separator = intersect( cluster.tree$clusters[[parent]],  cluster.tree$clusters[[child]]))
    cluster.tree$potentials[[parent]] <- abb[[2]]
    cluster.tree$potentials[[child]] <- abb[[1]]
  }

  clique.names <- names(V(cluster.tree$tree))

  cluster.tree$collected <- processed
  cluster.tree$active <- processed
  return(cluster.tree)
}

###########################################
## Distribute evidence
###########################################

DistributeEvidence <- function(cluster.tree, node){
  clique.names <- names(V(cluster.tree$tree))
  nodes_status <- data.frame(clique.names, active = FALSE, queued = FALSE, stringsAsFactors = FALSE) # are they actually inactive? NO, need to fix
  # but they actually are inactive, since propagate.worker clears cluster.tree$active!

  # nodes_status$active <- nodes_status$clique.names %in% cluster.tree$active
  nodes_status[nodes_status$clique.names %in% node,]$queued <- TRUE

  while(any(nodes_status$queued)){
    n <- nodes_status[nodes_status$queued,][1,1]
    nodes_status[nodes_status$queued,][1,]$active <- TRUE
    nodes_status[nodes_status$queued,][1,]$queued <- FALSE

    neighbors_names <- neighbors(cluster.tree$tree, n, mode = "all")$name # get neighbors of n
    neighbors_status <- nodes_status[nodes_status$clique.names %in% neighbors_names,] # get neighbors' rows
    inactive <- neighbors_status[!neighbors_status$active,] # cull for only uncollected neighbors
    inactive_names <- inactive$clique.names # get the names of uncollected neighbors

    if (length(inactive_names)>0) {
      nodes_status[nodes_status$clique.names %in% inactive$clique.names, ]$queued <- TRUE # queue inactive neighbors

      for (i in 1:length(inactive_names)) {
        abb <- Absorb(cluster.tree$potentials[[inactive_names[i]]], cluster.tree$potentials[[n]],
                      separator = intersect( cluster.tree$clusters[[n]],  cluster.tree$clusters[[inactive_names[i]]]),
                      distribute = TRUE)
        cluster.tree$potentials[[inactive_names[i]]] <- abb[[2]]
      }
    }

    cluster.tree$joint[[n]] <- cluster.tree$potentials[[n]]

  }
  cluster.tree$active <- nodes_status[nodes_status$active,]$clique.names

  return(cluster.tree)
}
