#' Plot the marginal distributions
#'
#' Plot the marginal distributions.
#'
#' @details Plot the marginal distributions. Marginals of discrete variables are plotted as
#' bar plots, while those of continuous variables as density plots.
#'
#' @param marginals the marginal distributions returned by \code{Marginals} for plotting
#' @param groups names of the marginals to be shown on plots
#'
#' @author Han Yu
#'
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks.
#' Journal of Machine Learning Research, 6(Sep), 1517-1550.
#'
#' @importFrom graphics lines legend barplot plot.default
#' @examples
#' data(toytree)
#' marg <- Marginals(toytree, c("Neu1", "Nr1i3", "chr1_42.65", "Spgl1"))
#' PlotMarginals(marginals=marg, groups=NULL)
#'
#' @seealso \code{\link{Marginals}}
#'
#' @export

PlotMarginals <- function(marginals, groups=NULL) {


  nms <- names(marginals$marginals)
  discrete.nodes <- nms[marginals$types]
  continuous.nodes <- nms[!marginals$types]

  posteriors <- marginals$marginals

  group.disc <- NULL
  group.cont <- NULL

  if(!is.null(groups)) {
    if(length(groups)!=length(posteriors)) {
      warning("Group and marginal lengths do not match.")
      groups <- NULL
    } else {
      group.disc <- groups[which(marginals$types)]
      group.cont <- groups[which(!marginals$types)]
    }
  }

  if(length(discrete.nodes)==0){
    PlotPosteriorContinuous(posteriors, groups=group.cont)
  }

  if(length(continuous.nodes)==0){
    par(mfrow=c(1,length(discrete.nodes)))
    for (i in 1:length(discrete.nodes)) {
      this.node <- discrete.nodes[i]
      PlotPosteriorDiscrete(posteriors[i], group=group.disc[i])
    }
    par(mfrow=c(1,1))
  }

  if(length(discrete.nodes)!=0 & length(continuous.nodes)!=0){
    par(mfrow=c(1,length(discrete.nodes)+1))
    PlotPosteriorContinuous(posteriors[continuous.nodes], groups=group.cont)
    for (i in 1:length(discrete.nodes)) {
      this.node <- discrete.nodes[i]
      PlotPosteriorDiscrete(posteriors[this.node], group=group.disc[i])
    }

    par(mfrow=c(1,1))
  }

}



