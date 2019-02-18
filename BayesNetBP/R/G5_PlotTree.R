#' Plot the cluster tree
#' 
#' Plot the structure of a \code{\linkS4class{ClusterTree}} object
#'
#' @details Plot the structure of \code{clustertree} object, with the nodes labeled by corresponding
#' elimination node. The circles represent continuous clusters, while the boxes represent discrete clusters.
#'
#' @param tree a \code{\linkS4class{ClusterTree}} object
#' @param color nodes color
#' 
#' @author Han Yu
#' @examples 
#' data(toytree)
#' PlotTree(toytree)
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550.  
#' 
#' @importFrom Rgraphviz plot
#' 
#' @export

PlotTree <- function(tree, color="gray90") {
  tree.graph <- tree@graph$tree
  cs.names <- tree@cluster
  cluster.type <- tree@cluster.class
  nAttrs <- list()
  node.shape <- c()
  node.shape[cluster.type] <- "box"
  node.shape[!cluster.type] <- "circle"
  names(node.shape) <- cs.names
  nAttrs$shape <- node.shape
  nAttrs$fontsize <- rep(16, length(cs.names))
  # nAttrs$height <- rep(1, length(cs.names))
  # nAttrs$width <- rep(2, length(cs.names))
  names(nAttrs$fontsize) <- cs.names
  # names(nAttrs$height) <- cs.names
  # names(nAttrs$width) <- cs.names
  fill.color <- rep(color, length(cs.names))
  names(fill.color) <- cs.names
  nAttrs$fillcolor <- fill.color
  Rgraphviz::plot(tree.graph, nodeAttrs=nAttrs, main="")
  
}

