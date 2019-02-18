#' Model compilation
#'
#' Compile the local models
#'
#' @details This function compiles the local models, including the conditional
#' probability tables for discrete variables, and linear predictor potentials
#' for continuous variables.
#' 
#' 
#' @param data a \code{data.frame} object or a \code{qtlnet} object
#' @param dag \code{NULL} if data is \code{qtlnet} object, or a \code{graphNEL} object of conditional 
#' Gaussian Bayesian network if data is \code{data.frame}.
#' @param node.class \code{NULL} if data is \code{qtlnet} object, or a \code{vector} of logical values 
#' named by node names, \code{TRUE} for discrete, \code{FALSE} for continuous variables if data 
#' is \code{data.frame}.
#' 
#' @return 
#' \describe{
#' \item{\code{pots}}{a \code{list} of discrete potentials (conditional probability tables) 
#' for each discrete variable. }
#' \item{\code{bags}}{a \code{list} of sets of continuous potentials (lppotentials), each set for a 
#' continuous variables.}
#' }
#'
#' @import doBy qtlnet
#' @importFrom qtl find.marker
#' @importFrom qtl scanone
#' @importFrom graph nodes
#' @importFrom igraph igraph.options
#' @author Han Yu
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. 
#' 
#' @examples 
#' 
#' data(liver)
#' models <- LocalModelCompile(data=liver$data, dag=liver$dag, node.class=liver$node.class)
#' 
#' @seealso \code{\link{ElimTreeInitialize}}
#'
#' @export

LocalModelCompile <- function(data, dag=NULL, node.class=NULL) {
  
  if("qtlnet" %in% class(data)){
    extr <- extractQTL(data)
    data <- extr$data
    dag <- extr$dag
    node.class <- extr$node.class
  }
  
  models <- ModelCompileData(data, dag, node.class)
  return(models)
}








