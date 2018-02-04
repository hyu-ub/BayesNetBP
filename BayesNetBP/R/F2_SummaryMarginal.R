#' Summary a continuous marginal distribution
#' 
#' This function summary the marginal distributions of continuous variables by outputing the
#' mean, standard deviation, and number of subpopulations
#' 
#' @param marginals the marginal distributions obtained from \code{\link{Marginals}} function
#' 
#' @return a \code{data.frame} object containing information about the marginal distributions for continuous variables. 
#' The marginal distributions of continous variables in a CG-BN model are mixtures of Gaussian distributions. 
#' Therefore, besides the mean and standard deviation, the object has an additional column to specify the number of Gaussian
#' mixtures.
#' 
#' \describe{
#'  \item{\code{mean}}{the mean value of a Gaussian distribution.}
#'  \item{\code{sd}}{the standard deviation of a Gaussian distribution.}
#'  \item{\code{n}}{the number of Gaussian distributions in the mixture.}
#' }
#' 
#' @examples 
#' 
#' data(liver)
#' tree.init.p <- Initializer(dag=liver$dag, data=liver$data, 
#'                            node.class=liver$node.class, 
#'                            propagate = TRUE)
#' marg <- Marginals(tree.init.p, c("HDL", "Ppap2a", "Neu1"))
#' SummaryMarginals(marginals=marg)
#' 
#' @seealso \code{\link{Marginals}}
#' 
#' @export

SummaryMarginals <- function(marginals) {
  
  margs <- marginals$marginals
  nms0 <- names(margs)
  types <- marginals$types
  
  nms <- mu <- sd <- n <-c()
  
  for (i in 1:length(types)) {
    if(!types[i]) {
      msd <- MeanSD(margs[[i]])
      nms <- c(nms, nms0[i])
      mu <- c(mu, msd[1])
      sd <- c(sd, msd[2])
      n <- c(n, nrow(margs[[i]]))
    }
  }
  
  df <- data.frame(Mean=mu, SD=sd, n=n)
  rownames(df) <- nms
  return(df)
}
