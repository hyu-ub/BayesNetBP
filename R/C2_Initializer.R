#' Initialize a ClusterTree object
#'
#' Initialize a ClusterTree object
#'
#' @details This is a wrapper function to initialize a \code{\linkS4class{ClusterTree}} object.
#'
#' @param dag a \code{graphNEL} object of the Bayesian network
#' @param data a \code{data.frame} object
#' @param node.class a named \code{vector} of \code{logical} values, \code{TRUE} if node 
#' is discrete, \code{FASLE} if otherwise
<<<<<<< HEAD
#' @param propagate \code{logical} \code{TRUE} if the discrete part of the \code{\linkS4class{ClusterTree}}
#' to be propagated
=======
>>>>>>> 391708a0cb06a04af75029f6c9bcc10c6c258fd4
#' 
#' @return \code{\linkS4class{ClusterTree}} object
#' 
#' @author Han Yu
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. 
#' 
#' @examples 
#' 
#' data(liver)
#' tree.init.p <- Initializer(dag=liver$dag, data=liver$data, node.class=liver$node.class)
#' @seealso The functions \code{\link{ClusterTreeCompile}} and \code{\link{LocalModelCompile}} provide necessary
#' 
#' objects to obtain \code{\linkS4class{ClusterTree}} object by initializing the elimination tree through this function.
#' 
#' @export

Initializer <- function(dag, data, node.class, propagate = TRUE){
  cst <- ClusterTreeCompile(dag=dag, node.class=node.class)
  models <- LocalModelCompile(data=data, dag=dag, node.class=node.class)
  tree.init <- ElimTreeInitialize(tree=cst$tree.graph, 
                                  dag=cst$dag, 
                                  model=models, 
                                  node.sets=cst$cluster.sets, 
                                  node.class=cst$node.class)
  if(propagate) {
<<<<<<< HEAD
    tree.init <- Propagate(tree.init)
=======
    tree.init <- PropagateDBN(tree.init)
>>>>>>> 391708a0cb06a04af75029f6c9bcc10c6c258fd4
  }
  return(tree.init)
}
