#' Absorb evidence into the model
#'
#' @details Absorb multiple types and pieces of evidences into a \code{\linkS4class{ClusterTree}}
#' object. The discrete compartment of the \code{\linkS4class{ClusterTree}} will be automatically
#' propagated after evidence absorption, so that the object will be ready for making
#' queries and absorbing additional evidence.
#'
#' @param tree a \code{\linkS4class{ClusterTree}} object
#' @param vars a \code{vector} of the names of observed variables
#' @param values a \code{list} of observed values of the variables. Aside from a single value,
#' The element of the list can also be a vector of likelihood values
#' 
#' @return \code{\linkS4class{ClusterTree}} object with the evidence absorbed
#' 
#' @author Han Yu
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. \cr
#' \cr
#' Lauritzen, S. L., & Spiegelhalter, D. J. (1988). Local computations with probabilities on 
#' graphical structures and their application to expert systems. Journal of the Royal Statistical 
#' Society. Series B (Methodological), 157-224.
#' 
#' @import stats graphics utils
#' @importFrom igraph igraph.from.graphNEL igraph.to.graphNEL V
#' 
#' 
#' @examples 
#' 
#' data(liver)
#' tree.init.p <- Initializer(dag=liver$dag, data=liver$data, 
#'                            node.class=liver$node.class, 
#'                            propagate = TRUE)
#' tree.post <- AbsorbEvidence(tree.init.p, c("Nr1i3", "chr1_42.65"), list(1,"1"))
#' 
#' @export

AbsorbEvidence <- function(tree, vars, values) {
  node.class <- tree@node.class
  
  hard <- c()
  soft <- c()
  hard.values <- list()
  soft.values <- list()
  
  if(sum(vars %in% tree@absorbed.variables)!=0) {
    var.in <- vars[vars %in% tree@absorbed.variables]
    msg1 <- paste0(var.in, collapse=", ")
    stop(paste0(msg1, " is/are already observed."))
  }
  
  if(sum(vars %in% tree@absorbed.soft.variables)!=0) {
    var.in <- vars[vars %in% tree@absorbed.soft.variables]
    msg1 <- paste0(var.in, collapse=", ")
    warning(paste0(msg1, " has/have absorbed likelihood evidence multiple times."))
  }
  
  if(length(vars)!=0){
    
    var.class <- node.class[vars]
    
    for(i in 1:length(vars)) {
      if (var.class[i]) {
        if (length(values[[i]])==1){
          tree <- DiscreteEvidence(tree, vars[i], values[[i]])
          hard <- c(hard, vars[i]) #
          hard.values <- append(hard.values, values[[i]]) #
        } 
        if (length(values[[i]])>1) {
          tree <- VirtualEvidence(tree, vars[i], values[[i]])
          soft <- c(soft, vars[i]) #
          soft.values <- append(soft.values, values[i]) #
        }
        
      } 
    }
    
    for(i in 1:length(vars)) {
      if (!var.class[i]) {
        tree <- PushOperation(tree, vars[i], values[[i]])
        hard <- c(hard, vars[i]) #
        hard.values <- append(hard.values, values[[i]]) #
      }
    }
    
  }
  
  tree <- Propagate(tree)
  
  tree@absorbed.variables <- c(tree@absorbed.variables, hard)
  tree@absorbed.soft.variables <- c(tree@absorbed.soft.variables, soft)
  tree@absorbed.values <- append(tree@absorbed.values, hard.values)
  tree@absorbed.soft.values <- append(tree@absorbed.soft.values, soft.values)
  names(tree@absorbed.values) <- tree@absorbed.variables
  names(tree@absorbed.soft.values) <- tree@absorbed.soft.variables
  
  return(tree)
}

