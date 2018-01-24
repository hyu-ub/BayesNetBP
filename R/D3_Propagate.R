
#' Propagate the cluster tree
#' 
#' This function propagates the discrete compartment of a \code{\linkS4class{ClusterTree}} object.
#' 
#' @details The discrete compartment must be propagted to get the joint distributions
#' of discrete variables in each discrete clusters. A \code{\linkS4class{ClusterTree}} object must be propagated
#' before absorbing evidence and making queries. 
#'
#' @param tree an initialized \code{\linkS4class{ClusterTree}} object
#' 
#' @return a \code{\linkS4class{ClusterTree}} object
#' 
#' @examples 
#' 
#' data(liver)
#' cst <- ClusterTreeCompile(dag=liver$dag, node.class=liver$node.class)
#' models <- LocalModelCompile(data=liver$data, dag=liver$dag, node.class=liver$node.class)
#' tree.init <- ElimTreeInitialize(tree=cst$tree.graph, 
#'                                 dag=cst$dag, 
#'                                 model=models, 
#'                                 node.sets=cst$cluster.sets, 
#'                                 node.class=cst$node.class)
#' tree.init@propagated
#' tree.init.p <- PropagateDBN(tree.init)
#' tree.init.p@propagated
#'
#' @export

Propagate <- function(tree) {
  
  discrete.clusters <- tree@cluster[tree@cluster.class]
  
  if (length(discrete.clusters)==0) {
    tree@propagated <- TRUE
    return(tree)
  }
  
  if (length(discrete.clusters)==1) {
    tree@jpt <- tree@cpt
    tree@propagated <- TRUE
    return(tree)
  }
  
  tree.graph <- igraph.from.graphNEL(tree@graph$tree)
  tree.sub.graph <- induced_subgraph(tree.graph, discrete.clusters)
  potentials.sub <- tree@cpt[discrete.clusters]
  discrete.sets <- tree@member[discrete.clusters]
  tree@jpt <- propagate.worker(tree.sub.graph, potentials.sub, discrete.sets)
  tree@propagated <- TRUE
  return(tree)
  
}
