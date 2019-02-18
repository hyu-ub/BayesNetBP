
#' Launch the BayesNetBP Shiny App
#'
#' @details The function \code{runBayesNetApp} lauches the \code{Shiny App} accompanied
#' with this package. The app loads the \code{toytree} example by default and allows users
#' to load customized \code{\linkS4class{ClusterTree}} object. In order to use this feature, a 
#' \code{\linkS4class{ClusterTree}} object should be built, propagated and named \code{tree.init.p}, and 
#' then saved as a \code{.RDATA} file. This file can be read in by the app. \cr
#' 
#' The console of \code{BayesNetBP} Shiny App comprises three panels. The first 
#' part controls the model loading and network layouts.  It also allows user to subset the 
#' network to faciliate visualization. The \code{Expand} function can trace the ancestors, descendants, or both, of 
#' a selected node in a stepwise manner. The expanded nodes will be colored orange. By clicking 
#' \code{Add to list}, the expanded nodes will be selected and be purple. The user can continue
#' selecting other nodes by using \code{Expand} and \code{Add to list} functions at this stage. After
#' selecting desired node sets, the user can subset the graph by the \code{Subset} function. The nodes in
#' subsetted graph retain all properties before subsetting, including their colors and divergence. \cr
#' 
#' The second panel is used for absorption of fixed and hard evidences.  
#' The users can add multiple pieces of evidence to a list and absorb them into the model simultaneously.  
#' The nodes with evidence absorbed will be colored green when the absorption is complete. Marginals of the nodes can be quried as 
#' density or bar plots by node types. If a set of evidence has been absorbed, the marginals both 
#' before and after absorption will be returned to facilitate comparison. To query the marginals, the user can 
#' select the node of interest in the graph, and then click \code{Plot Marginals}. The \code{Shift in Marginals} 
#' function computes the signed and symmetric Kullback-Liebler divergence for all applicable nodes 
#' in the network, and colors the nodes in a similar manner as the function \code{\link{PlotCGBN}}. \cr
#' 
#' The function for systematic assessment of variable marginal shifts is provided in the third panel. 
#' It allows user to specify which node to absorb the spectrum of evidence in a menu, and to select whose 
#' divergence to be calculated by firstly selecting the node on the graph and then clicking \code{Add to Plot List}.
#' Alternatively, the user can use \code{Add All} function to select all applicable nodes into the plotting list.
#' The result is visualized in an interactive plot.
#'
#' @param launch.browser \code{logical(1)} whether launch the App in browser
#' 
#' @author Han Yu
#' 
#' @examples 
#' 
#' \dontrun{
#' # load or install required packages to run App
#' library("shiny")
#' library("googleVis")
#' library("devtools")
#' devtools::install_github("cytoscape/r-cytoscape.js")
#' # run the App in browser
#' runBayesNetApp(launch.browser=TRUE)
#' }
#' 
#' @importFrom igraph igraph.from.graphNEL as_edgelist
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