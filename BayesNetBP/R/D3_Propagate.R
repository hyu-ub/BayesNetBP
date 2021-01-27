
#' Propagate the cluster tree
#'
#' This function propagates the discrete compartment of a \code{\linkS4class{ClusterTree}} object.
#'
#' @details The discrete compartment must be propagted to get the joint distributions
#' of discrete variables in each discrete clusters. A \code{\linkS4class{ClusterTree}} object must be propagated
#' before absorbing evidence and making queries.
#'
#' @param tree an initialized \code{\linkS4class{ClusterTree}} object
#' @param targets the cluster involved in evidence propagation, usually set by default
#'
#' @return a \code{\linkS4class{ClusterTree}} object
#'
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks.
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. \cr
#' \cr
#' Lauritzen, S. L., & Spiegelhalter, D. J. (1988). Local computations with probabilities on
#' graphical structures and their application to expert systems. Journal of the Royal Statistical
#' Society. Series B (Methodological), 157-224. \cr
#' \cr
#' Yu H, Moharil J, Blair RH (2020). BayesNetBP: An R Package for Probabilistic Reasoning in Bayesian
#' Networks. Journal of Statistical Software, 94(3), 1-31. <doi:10.18637/jss.v094.i03>.
#'
#' @examples
#'
#' data(liver)
#' tree.init <- Initializer(dag=liver$dag, data=liver$data,
#'                          node.class=liver$node.class,
#'                          propagate = FALSE)
#' tree.init@propagated
#' tree.init.p <- Propagate(tree.init)
#' tree.init.p@propagated
#'
#' @export

# targets is a vector of any of the nodes in the graph, or NA
Propagate <- function(tree, targets = NA) {

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
  tree@jpt <- propagate.worker(tree.sub.graph, potentials.sub, discrete.sets, targets = targets)
  tree@propagated <- TRUE
  return(tree)

}
