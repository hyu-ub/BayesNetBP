
#' Launch the BayesNetBP Shiny App
#'
#' @details The function \code{runBayesNetApp} lauches the \code{Shiny App} accompanied
#' with this package. The app loads the \code{toytree} example by default and allows users
#' to load customized \code{\linkS4class{ClusterTree}} object. In order to use this feature, a
#' \code{\linkS4class{ClusterTree}} object should be built, propagated and named \code{tree.init.p}, and
#' then saved as a \code{.RDATA} file. This file can be read in by the app. \cr
#'
#' The console of \code{BayesNetBP} Shiny App comprises three panels. The first
#' part controls the model loading, visualization and subnetwork selection. The \code{Fit} function fits
#' the entire graph in the window. The \code{Fit Selected} function fits the selected subnetwork to the window.
#' The user can subset the network for visualization. The \code{Expand} function can trace the one hop neighbor of
#' selected nodes in a stepwise manner.
#' After selecting desired node sets, the user can subset the graph by the \code{Subset} function.
#' \cr
#'
#' The second panel is used for absorption of fixed and hard evidences.
#' The users can add multiple pieces of evidence to a list and absorb them into the model simultaneously.
#' Marginals of other nodes can be quried as
#' density or bar plots by node types. If a set of evidence has been absorbed, the marginals both
#' before and after absorption will be returned to facilitate comparison. To query the marginals, the user can
#' select the node of interest in the graph, and then click \code{Marginal of Selected}. The \code{Shift in Marginals}
#' function computes the signed and symmetric Kullback-Liebler divergence for all applicable nodes
#' in the network, and colors the nodes by their divergence and change in directions. \cr
#'
#' The function for systematic assessment of variable marginal shifts is provided in the third panel.
#' It allows user to specify which node to absorb the spectrum of evidence in the select menu and click \code{Select Observed}, and to select whose
#' divergence to be calculated by selecting the node in the menu and then clicking \code{Add to Plot}.
#' Alternatively, the user can use \code{Add All} function to select all applicable nodes into the plotting list.
#' The result is visualized in an interactive plot. The \code{Min}, \code{Max} and \code{Step} controls the range of values
#' of the evidence to be absorbed.
#'
#' @param launch.browser \code{logical(1)} whether launch the App in browser
#'
#' @author Han Yu
#'
#' @references Yu H, Moharil J, Blair RH (2020). BayesNetBP: An R Package for Probabilistic Reasoning in Bayesian
#' Networks. Journal of Statistical Software, 94(3), 1-31. <doi:10.18637/jss.v094.i03>.
#'
#' @examples
#'
#' \dontrun{
#' # load or install required packages to run App
#' library("shiny")
#' library("googleVis")
#' library("devtools")
#' devtools::install_github("cytoscape/cyjShiny")
#' library("cyjShiny")
#' # run the App in browser
#' runBayesNetApp(launch.browser=TRUE)
#' }
#'
#' @importFrom igraph igraph.from.graphNEL as_edgelist
#'
#' @references Yu H, Moharil J, Blair RH (2020). BayesNetBP: An R Package for Probabilistic Reasoning in Bayesian
#' Networks. Journal of Statistical Software, 94(3), 1-31. <doi:10.18637/jss.v094.i03>.
#'
#' @export

runBayesNetApp <- function(launch.browser=TRUE) {
  appDir <- system.file("BayesNetApp", package = "BayesNetBP")
  if (appDir == "") {
    stop("Could not find BayesNetApp directory. Try re-installing `BayesNetBP`.", call. = FALSE)
  }
  data(toytree, envir = environment())
  shiny::runApp(appDir, launch.browser=launch.browser)
}
