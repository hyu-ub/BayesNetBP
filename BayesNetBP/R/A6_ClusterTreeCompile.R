#' Compile the cluster tree
#'
#' Get the cluster sets and strong semi-elimination tree from the Bayesian network
#'
#' @details This function forms the cluster sets and the semi-elimination tree graph
#' from the Bayesian network. The procedures include acquiring the elimination order,
#' moralization, triangulation, obtaining cluster sets, forming strong elimination 
#' tree and strong semi-elimination tree. The cluster sets and the semi-elimination 
#' tree are required to initialize the cluster tree.
#'
#' @param dag a \code{graphNEL} object of the Bayesian network
#' @param node.class a named \code{vector} of \code{logical} values, \code{TRUE} if node 
#' is discrete, \code{FASLE} if otherwise
#' @return 
#' \describe{
#' \item{\code{tree.graph}}{a \code{graphNEL} object of semi-elimination tree.}
#' \item{\code{dag}}{a \code{graphNEL} object of original Bayesian network.}
#' \item{\code{cluster.sets}}{a \code{list} of members of each cluster.}
#' \item{\code{node.class}}{a named \code{vector} of \code{logical} values, \code{TRUE} if node 
#' is discrete, \code{FASLE} if otherwise}
#' \item{\code{elimination.order}}{a \code{vector} of node names sorted by the elimination order.}
#' }
#' 
#' @author Han Yu
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. 
#'
#' @importFrom igraph neighbors add_edges are_adjacent delete_vertices graph_from_adjacency_matrix graph_from_data_frame topological.sort
#' 
#' @examples 
#' 
#' data(liver)
#' cst <- ClusterTreeCompile(dag=liver$dag, node.class=liver$node.class)
#' 
#' @seealso \code{\link{ElimTreeInitialize}}
#' 
#' @export

ClusterTreeCompile <- function(dag, node.class) {
  
  elim.order <- EliminationOrder(dag, node.class=node.class)
  graph.mor <- Moralize(dag)
  graph.tri <- Triangulate(graph.mor, elim.order)
  cs <- ElimTreeNodes(graph.tri, elim.order)
  strongET <- StrongEliminationTree(cs, elim.order)
  semiET <- SemiEliminationTree(strongET, cs, node.class, elim.order)
  output <- list(tree.graph=semiET,
                 dag=dag,
                 cluster.sets=cs,
                 node.class=node.class,
                 elimination.order=elim.order)
  return(output)
  
}



