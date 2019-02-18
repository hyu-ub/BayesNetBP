#' Possible values of a discrete variable
#'
#' Obtain all the possible values of a discrete variable.
#'
#' @param tree a \code{\linkS4class{ClusterTree}} object
#' @param var the variables to be queried
#' @param message type of desired distribution
#' @return a \code{vector} of the possible values of discrete variable. If the variable is continuous, 
#' the returned value will be \code{NULL}.
#' 
#' @examples 
#' data(toytree)
#' GetValue(toytree, "HDL")
#' 
#' @author Han Yu
#' 
#' @export

GetValue <- function(tree, var, message=TRUE) {
  
  if (!var %in% tree@node) {
    if (message) {
      cat("This node is not found. \n")
    }
    return(NULL)
  }
  
  if (!tree@node.class[var]) {
    if (message) {
      cat("The node is continuous. \n")
    }
    return(NULL)
  }
  
  if (var %in% tree@absorbed.variables) {
    if (message) {
      cat("The node is absorbed with value", tree@absorbed.values[[var]], "\n")
    }
    return(NULL)
  }
  
  if (var %in% tree@cluster) {
    tab <- tree@cpt[[var]]$cpt
    clname <- colnames(tab)
    var.ind <- which(clname==var)
    return(as.character(unique(tab[,var.ind])))
  }
  
  for (i in 1:length(tree@cpt)) {
    tab <- tree@cpt[[i]]$cpt
    clname <- colnames(tab)
    if (var %in% clname) {
      var.ind <- which(clname==var)
      return(as.character(unique(tab[,var.ind])))
    }
  }
  
}