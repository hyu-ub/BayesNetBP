#' Initialize the elimination tree
#'
#' Initialize the elimination tree with the local models
#'
#' @details Initialize the elimination tree with the local models
#'
#' @param tree a \code{graphNEL} object of the elimination tree
#' @param dag a \code{graphNEL} object of the Bayesian network
#' @param model a \code{list} of local models built from \code{\link{LocalModelCompile}} function
#' @param node.sets a \code{list} of cluster sets obtained from \code{\link{ClusterTreeCompile}} function
#' @param node.class a named \code{vector} of \code{logical} values, \code{TRUE} if node
#' is discrete, \code{FASLE} if otherwise
#'
#' @return \code{\linkS4class{ClusterTree}} object with the local models incorporated
#'
#' @author Han Yu
#'
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks.
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. \cr
#' \cr
#' Yu H, Moharil J, Blair RH (2020). BayesNetBP: An R Package for Probabilistic Reasoning in Bayesian
#' Networks. Journal of Statistical Software, 94(3), 1-31. <doi:10.18637/jss.v094.i03>.
#'
#' @import doBy
#' @importFrom graph nodes
#' @importFrom igraph neighbors
#' @importFrom methods new
#' @examples
#'
#' data(liver)
#' cst <- ClusterTreeCompile(dag=liver$dag, node.class=liver$node.class)
#' models <- LocalModelCompile(data=liver$data, dag=liver$dag, node.class=liver$node.class)
#' tree.init <- ElimTreeInitialize(tree=cst$tree.graph,
#'                                 dag=cst$dag,
#'                                 model=models,
#'                                 node.sets=cst$cluster.sets,
#'                                 node.class=cst$node.class)
#'
#' @seealso The functions \code{\link{ClusterTreeCompile}} and \code{\link{LocalModelCompile}} provide necessary
#' objects to obtain \code{\linkS4class{ClusterTree}} object by initializing the elimination tree through this function.
#'
#' @export

ElimTreeInitialize <- function(tree, dag, model, node.sets, node.class){

  e.seq <- EliminationOrder(dag, node.class)
  ClusterTree <- new("ClusterTree",
                     cluster = graph::nodes(tree),
                     node = e.seq,
                     graph = list(dag = dag, tree = tree),
                     member = node.sets,
                     node.class = node.class[e.seq],
                     propagated = FALSE
                     )

  ClusterTree@activeflag <- rep(TRUE, length(ClusterTree@cluster))
  names(ClusterTree@activeflag) <- ClusterTree@cluster

  tree.graph <- igraph.from.graphNEL(tree) # get the igraph object for tree
  dag.graph <- igraph.from.graphNEL(dag)

  for (i in 1:length(ClusterTree@cluster)) {

    this.cluster <- ClusterTree@cluster[i]
    this.par <- neighbors(tree.graph, v=this.cluster, mode="in")$name
    if (length(this.par)==0) {
      ClusterTree@parent[i] <- NA
    } else {
      ClusterTree@parent[i] <- this.par
    }

    ClusterTree@cluster.class[i] <- as.logical(prod(node.class[ node.sets[[ ClusterTree@cluster[i] ]] ]))

  }

  names(ClusterTree@parent) <- ClusterTree@cluster
  names(ClusterTree@cluster.class) <- ClusterTree@cluster

  continuous.clusters <- ClusterTree@cluster[!ClusterTree@cluster.class]
  continuous.nodes <- names(node.class)[!node.class]
  discrete.clusters <- ClusterTree@cluster[ClusterTree@cluster.class]
  discrete.nodes <- names(node.class)[node.class]

  ClusterTree@assignment <- asgn <- assignUniverse(dag, node.sets[discrete.clusters], discrete.nodes)

  #################################################
  ## initialize with local models
  #################################################

  ## initialize the discrete part

  if(length(discrete.clusters)!=0) {

    for (i in 1:length(discrete.clusters)) {
      this.cluster <- discrete.clusters[i]
      for (j in 1:length(asgn[[this.cluster]])) {
        this.asgn <- asgn[[this.cluster]][j]
        if (j==1) {
          pot <- model$pots[[this.asgn]]
        } else {
          pot <- factor.product(pot, model$pots[[this.asgn]])
        }
      }
      ClusterTree@cpt[[i]] <- pot
    }
    names(ClusterTree@cpt) <- discrete.clusters

  }

  ## initialize the continuous part

  if (length(continuous.clusters)!=0) {

    for (j in 1:length(continuous.clusters)) {
      ClusterTree@lppotential[[j]] <- list()
      ClusterTree@postbag[[j]] <- list()
    }

    names(ClusterTree@lppotential) <- continuous.clusters
    names(ClusterTree@postbag) <- continuous.clusters

    reallocate <- FALSE
    cl.reallocate <- c()

    for (i in 1:length(continuous.nodes)) {
      this.node <- continuous.nodes[i]
      this.par <- neighbors(dag.graph, v=this.node, mode="in")$name
      this.all <- c(this.node, this.par)
      for (j in 1:length(continuous.clusters)) {
        this.cluster <- continuous.clusters[j]
        this.member <- ClusterTree@member[[this.cluster]]
        if (is.subset(this.all, this.member)) {
          if(this.node==this.cluster){
            l <- length(ClusterTree@lppotential[[j]])
            ClusterTree@lppotential[[j]][[l+1]] <- model$bags[[this.node]]
            names(ClusterTree@lppotential[[j]]) <- c(names(ClusterTree@lppotential[[j]]), this.node)
          } else {
            reallocate <-TRUE # need reallocation of lppotentials
            cl.reallocate <- union(cl.reallocate, this.cluster) # record the clusters requiring reallocation
            l <- length(ClusterTree@postbag[[j]])
            ClusterTree@postbag[[j]][[l+1]] <- model$bags[[this.node]]
            names(ClusterTree@postbag[[j]]) <- c(names(ClusterTree@postbag[[j]]), this.node)
          }
          break
        }
      }
    }

    # This step seems to be unnecessary if the elimination sequence is properly arranged.
    # No reallocation of lppotentials in postbags is needed in the data sets tested so far.
    # In case of finding such a case, please report it to the maintainer

    if (reallocate) {
      cat("Reallocation of LPPotentials is required for clusters: ",
          paste0(cl.reallocate, collapse=", "), ".", "Please contact maintainer.\n")
      # Reallocate
    }

  }
  return(ClusterTree)
}


