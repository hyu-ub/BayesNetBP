###########################################
## Propagate
###########################################

propagate.worker <- function(tree.graph, potentials, cluster.sets){
  
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

CollectEvidence <- function(cluster.tree, node) {
  clique.names <- names(V(cluster.tree$tree))
  ngbs <- neighbors(cluster.tree$tree, node, mode = "all")$name
  inactive <- setdiff(ngbs, cluster.tree$active)
  cluster.tree$active <- c(cluster.tree$active, node)
  
  for (ngb in inactive){
    # cat("collecting for ", node, "from", ngb)
    cluster.tree <- CollectEvidence(cluster.tree, ngb)
    # collected <- ce[[1]]
    # cst <- ce[[2]]
  }
  
  if (length(inactive)>0) {
    for (i in 1:length(inactive)) {
      abb <- Absorb(cluster.tree$potentials[[node]], cluster.tree$potentials[[inactive[i]]],
                    separator = intersect( cluster.tree$clusters[[node]],  cluster.tree$clusters[[inactive[i]]]))
      cluster.tree$potentials[[node]] <- abb[[2]]
      cluster.tree$potentials[[inactive[i]]] <- abb[[1]]
      # cat(inactive[i], " -> ", node, "\n")
    }
    
  }
  cluster.tree$collected <- c(cluster.tree$collected, node)
  return(cluster.tree)
}

###########################################
## Distribute evidence
###########################################

# need to reset the active nodes of cluster.tree after collecting evidence

DistributeEvidence <- function(cluster.tree, node){
  clique.names <- names(V(cluster.tree$tree))
  ngbs <- neighbors(cluster.tree$tree, node, mode = "all")$name
  
  inactive <- setdiff(ngbs, cluster.tree$active)
  cluster.tree$active <- c(cluster.tree$active, node)
  
  cluster.tree$joint[[node]] <- cluster.tree$potentials[[node]]
  
  if (length(inactive)>0) {
    for (i in 1:length(inactive)) {
      abb <- Absorb(cluster.tree$potentials[[inactive[i]]], cluster.tree$potentials[[node]],
                    separator = intersect( cluster.tree$clusters[[node]],  cluster.tree$clusters[[inactive[i]]]),
                    distribute = TRUE)
      cluster.tree$potentials[[inactive[i]]] <- abb[[2]]
      cluster.tree <- DistributeEvidence(cluster.tree, inactive[i])
    }
  }
  
  return(cluster.tree)
}

