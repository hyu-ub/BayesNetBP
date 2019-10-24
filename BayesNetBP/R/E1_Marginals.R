#' Obtain marginal distributions
#'
#' Get the marginal distributions of multiple variables
#'
#' @details Get the marginal distributions of multiple variables. The function \code{Marginals}
#' returns a \code{list} of marginal distributions. The marginal distribution of a discrete variable
#' is a named vector of probabilities. Meanwhile, the marginal distributions of
#' continous variables in a CG-BN model are mixtures of Gaussian distributions.
#' To fully represent this information, the marginal of a continuous variable is represented by
#' a \code{data.frame} with three columns to specify
#' parameters for each Gaussian distribution in the mixture, which are
#'
#' \describe{
#'  \item{\code{mean}}{the mean value of a Gaussian distribution.}
#'  \item{\code{sd}}{the standard deviation of a Gaussian distribution.}
#'  \item{\code{n}}{the number of Gaussian mixtures}
#' }
#'
#' @param tree a \code{\linkS4class{ClusterTree}} object
#' @param vars a \code{vector} of variables for query of marginal distributions
#'
#' @return
#'
#' \describe{
#'  \item{\code{marginals}}{a \code{list} of marginal distributions}
#'  \item{\code{types}}{a named \code{vector} indicating the types of the variables whose
#'  marginals are queried: \code{TRUE} for discrete, \code{FALSE} for continuous.}
#' }
#'
#' @author Han Yu
#'
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks.
#' Journal of Machine Learning Research, 6(Sep), 1517-1550.
#'
#' @examples
#'
#' data(liver)
#' tree.init.p <- Initializer(dag=liver$dag, data=liver$data,
#'                            node.class=liver$node.class,
#'                            propagate = TRUE)
#' tree.post <- AbsorbEvidence(tree.init.p, c("Nr1i3", "chr1_42.65"), list(1,"1"))
#' marg <- Marginals(tree.post, c("HDL", "Ppap2a"))
#' marg$marginals$HDL
#' head(marg$marginals$Ppap2a)
#'
#' @seealso \code{\link{PlotMarginals}} for visualization of the marginal distributions,
#' \code{\link{SummaryMarginals}} for summarization of the marginal distributions of
#' continuous variables.
#'
#' @export

Marginals <- function(tree, vars) {

  if(!tree@propagated) {
    stop("The ClusterTree object must be propagated before making queries.")
  }

  if(sum(vars %in% tree@absorbed.variables)!=0) {
    var.in <- vars[vars %in% tree@absorbed.variables]
    msg1 <- paste0(var.in, collapse=", ")
    stop(paste0(msg1, " is/are already observed."))
  }

  node.class <- tree@node.class
  marginal.types <- node.class[vars]

  margs <- list()

  for (i in 1:length(vars)) {
    var <- vars[i]
    if (node.class[[var]]) {
      margs[[i]] <- DiscreteMarginal(tree, var)
    } else {
      margs[[i]] <- PushMarginal(tree, var)
    }
  }

  names(margs) <- vars
  output <- list()
  output$marginals <- margs
  output$types <- marginal.types
  class(output) <- "marginals"

  return(output)
}

