#' An S4 class of the cluster tree.
#' 
#' @description The \code{ClusterTree} object is the computational object for belief propagation.
#' @slot cluster A \code{vector} storing the name of clusters in the cluster tree.
#' @slot node A \code{vector} storing the name of nodes in the Bayesian network.
#' @slot graph A \code{list} of two graphNEL objects: \code{$dag} stores the graph of Bayesian network,
#'               \code{$tree} stores the graph of the cluster tree.
#' @slot member A named \code{list} of the node cluster membership.
#' @slot parent A named \code{vector} indicating the parent node of a given cluster in the cluster tree.
#' @slot cluster.class A named \code{vector} of logical values indicating whether a cluster is continuous or discrete.
#' @slot node.class A named \code{vector} of logical values indicating whether a node is continuous or discrete. 
#' @slot assignment A named \code{list} indicating the assignment of discrete nodes discrete clusters.
#' @slot propagated A \code{logical} value indicating whether the discrete compartment has been propagated. 
#' 
#' @slot cpt A named \code{list} of the conditional probability tables.
#' @slot jpt A named \code{list} of the joint distribution tables.
#' @slot lppotential A named \code{list} of the linear predictor potentials assigned to each cluster in the lppotential slots.
#' @slot postbag A named \code{list} of the linear predictor potentials assigned to each cluster in the postbag slots.
#' @slot activeflag A named \code{vector} of logical values indicating whether a continuous cluster is active.

#' @slot absorbed.variables A \code{vector} of characters indicating variables observed with hard evidence.
#' @slot absorbed.values A \code{list} indicating the values of the variables observed with hard evidence. 
#' @slot absorbed.soft.variables A \code{vector} of characters indicating variables observed with soft or likelihood evidence.
#' @slot absorbed.soft.values A \code{list} of the likelihoods of the soft or likelihood evidence.

setClass("ClusterTree",
         slots = list(cluster = "character", 
                      node = "character",
                      graph = "list",
                      member = "list",
                      parent = "character",
                      cluster.class = "logical",
                      node.class = "logical",
                      assignment = "list",
                      propagated = "logical",
                      
                      cpt = "list",
                      jpt = "list",
                      lppotential = "list",
                      postbag = "list",
                      activeflag = "logical",
                      
                      absorbed.variables = "character",
                      absorbed.values = "list",
                      absorbed.soft.variables = "character",
                      absorbed.soft.values = "list"
         )
)

