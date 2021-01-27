#' Convert qtlnet to graphNEL object
#'
#' Extract network structure from qtlnet object and convert to graphNEL object
#'
#' @details This function extracts network structure from qtlnet object and convert to graphNEL object.
#' The example data can be downloaded from <https://github.com/hyu-ub/BayesNetBP>.
#' @param data a \code{qtlnet} object
#'
#' @return \item{\code{graphNEL}}{a \code{graphNEL} object.}
#'
#' @importFrom igraph igraph.to.graphNEL
#' @author Han Yu
#'
#' @examples
#'
#' \dontrun{
#' load(liverqtl.rda)
#' qtlnet_to_graphNEL(liverqtl$qtlnet.fit)
#' }
#'
#' @export

qtlnet_to_graphNEL <- function(data) {
  qtlnet_igraph <- igraph_from_qtlnet(data)
  qtlnet_graphnel <- igraph::igraph.to.graphNEL(qtlnet_igraph)
  return(qtlnet_graphnel)
}



