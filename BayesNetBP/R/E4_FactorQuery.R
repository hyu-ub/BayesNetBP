#' Queries of discrete variable distributions 
#'
#' Obtain the joint, marginal, and conditional distributions of discrete variables
#'
#' @details Query the joint distribution of any combination of discrete variables when 
#' mode is "joint", or conditional distribution of a discrete variable. The mode "list"
#' return a \code{list} of variable combinations, such that joint distributions of any subset 
#' of them are ready for extraction. Queries outside this list are also supported but may 
#' take longer computing time. This function will also return marginal distribution if only
#' one variable is queried.
#' 
#' @param tree a \code{\linkS4class{ClusterTree}} object
#' @param vars the variables to be queried
#' @param mode type of desired distribution
#' @return \code{data.frame} object specifying a joint or conditional distribution.
#' 
#' @author Han Yu
#' 
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks. 
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. 
#' 
#' @importFrom igraph neighbors all_simple_paths induced_subgraph
#' 
#' @examples 
#' 
#' data(chest)
#' dag <- chest$dag
#' node.class <- rep(TRUE, length(dag@nodes))
#' names(node.class) <- dag@nodes
#' tree.init.p <- Initializer(dag=dag, data=chest$data, 
#'                            node.class=node.class, 
#'                            propagate=TRUE)
#' # joint distribution
#' FactorQuery(tree=tree.init.p, vars=c("tub", "xray", "dysp", "asia"), mode="joint")
#' 
#' # conditional distribution
#' FactorQuery(tree=tree.init.p, vars=c("xray"), mode="conditional")
#' 
#' @export

FactorQuery <- function(tree, vars=c(), mode=c("joint", "conditional", "list")) {
  
  if (length(tree@jpt)==0) stop("Joint distribution table not found. 
                                No discrete nodes in model or the ClusterTree is not propagated.")
  
  if(sum(vars %in% tree@absorbed.variables)!=0) {
    var.in <- vars[vars %in% tree@absorbed.variables]
    msg1 <- paste0(var.in, collapse=", ")
    stop(paste0(msg1, " is/are already observed."))
  }
  
  if (mode=="list") {
    result <- list()
    j <- 1
    for (i in 1:length(tree@jpt)) {
      if(is.subset(vars, colnames(tree@jpt[[i]]$cpt))) {
        result[[j]] <- colnames(tree@jpt[[i]]$cpt)
        j <- j+1
      }
    }
    return(result)
  }
  
  dag.graph <- igraph.from.graphNEL(tree@graph$dag)
  
  if (mode=="conditional") {
    if (length(vars)!=1) {
      stop("If mode is conditional, vars should be a single variable.")
    }
    parents <- names(neighbors(dag.graph, vars, mode="in"))
    parents <- setdiff(parents, tree@absorbed.variables)
    allvar <- c(vars, parents)
    for (i in 1:length(tree@jpt)){
      if(is.subset(allvar, colnames(tree@jpt[[i]]$cpt))) {
        pot <- tree@jpt[[i]]
        result <- marginalize.discrete(pot, allvar)
        result <- conditional(result, parents)
        output <- data.frame(result$cpt, prob=result$prob)
        break
      }
    }
    rownames(output) <- NULL
    return(output)
  }
  
  
  if (mode=="joint") {
    for (i in 1:length(tree@jpt)) {
      if(is.subset(vars, colnames(tree@jpt[[i]]$cpt))) {
        pot <- tree@jpt[[i]]
        result <- marginalize.discrete(pot, vars)
        output <- data.frame(result$cpt, prob=result$prob)
        rownames(output) <- NULL
        return(output)
      }
    }
    result <- query.ooc(tree, vars)
    output <- data.frame(result$cpt, prob=result$prob)
    rownames(output) <- NULL
    return(output)
  }
  
}