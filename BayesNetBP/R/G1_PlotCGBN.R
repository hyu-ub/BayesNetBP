#' Plot the Bayesian network
#'
#' Plot and compare two Bayesian networks with different evidence(s) absorbed and propagated.
#'
#' @details Network visualization of the node-specific differences between Bayesian Networks
#' with the same topology, but evidence that has been absorbed and propagated.  The change of
#' marginal distribution of each node is measured by signed and symmetric Kullback-Leibler
#' divergence.  The sign indicates the direction of change, with \code{tree.1} considered as the baseline.
#' The magnitude of the change is reflected by the value.  Nodes that are white are d-separated
#' from the evidence. This function requires \code{Rgraphviz} package.
#'
#' @param tree.1 a \code{\linkS4class{ClusterTree}}
#' @param tree.2 a \code{\linkS4class{ClusterTree}}
#' @param fontsize font size for the node labels
#' @param pbar \code{logical(1)} whether to show progress bar
#' @param plotting \code{logical(1)} whether to output plot
#' @param epsilon \code{numeric(1)} the KL divergence is undefined if certain states of a discrete variable 
#' have probabilities of 0. In this case, a small positive number epsilon is assigned as their probabilities for calculating
#' the divergence. The probabilities of other states are shrunked proportionally to ensure they sum up to 1.
#' @return a plot of Bayesian network
#' @return a \code{vector} of signed symmetric Kullback-Leibler divergence
#'
#' @import RColorBrewer
#'
#' @author Han Yu
#'
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks.
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. \cr
#' \cr
#' Yu H, Moharil J, Blair RH (2020). BayesNetBP: An R Package for Probabilistic Reasoning in Bayesian
#' Networks. Journal of Statistical Software, 94(3), 1-31. <doi:10.18637/jss.v094.i03>.
#'
#' @importFrom fields image.plot
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics par
#'
#' @examples
#' \dontrun{
#' library("Rgraphviz")
#' data(toytree)
#' tree.post <- AbsorbEvidence(toytree, c("Nr1i3"), list(1))
#' PlotCGBN(tree.1=toytree, tree.2=tree.post)
#' }
#' @export

PlotCGBN <- function(tree.1, tree.2, fontsize=NULL, pbar=FALSE, plotting=TRUE, epsilon = 10^-6) {

  # tree.1 <- tree.init.p; tree.2 <- tree.post; fontsize=NULL; pbar=TRUE; plotting=TRUE;
  dag <- tree.1@graph$dag
  node.class <- tree.1@node.class

  absorbed.1 <- tree.1@absorbed.variables
  absorbed.2 <- tree.2@absorbed.variables
  absorbed <- union(absorbed.1, absorbed.2)

  var.inter <- intersect(absorbed.1, absorbed.2)
  var.1 <- setdiff(absorbed.1, var.inter)
  var.2 <- setdiff(absorbed.2, var.inter)

  vars <- absorbed
  node.names <- tree.1@node
  active.vars <- setdiff(node.names, vars)

  disc.all <- names(node.class)[node.class]
  cont.all <- names(node.class)[!node.class]

  disc.active <- intersect(disc.all, active.vars)
  cont.active <- intersect(cont.all, active.vars)

  all.active <- c(disc.active, cont.active)

  klds.0 <- c()

  sys <- Sys.info()[1]

  if(pbar){
    if(sys=="Windows") {
      pb <- winProgressBar(title = "Computing posterior", min = 0,
                           max = length(all.active), width = 300)
    } else {
      pb <- txtProgressBar(min = 0, max = length(all.active), style = 3)
    }
  }


  for (i in 1:length(all.active)){

    if(pbar){
      if(sys=="Windows") {
        setWinProgressBar(pb, i, title=paste("Computing posterior for ", all.active[i], ": ",
                                             round((i-1)/length(all.active)*100, 0), "% complete"))
      } else {
        setTxtProgressBar(pb, i)
      }
    }

    post.1 <- Marginals(tree.1, all.active[i])
    post.2 <- Marginals(tree.2, all.active[i])
    klds.0[i] <- SymmetricKLD(post.1[[1]][[1]], post.2[[1]][[1]], 
                              discrete = node.class[all.active[i]], epsilon = epsilon) ##
  }
  names(klds.0) <- all.active
  klds <- abs(klds.0)
  if(pbar){close(pb)}

  if (!plotting){
    return(klds.0)
  }

  nAttrs <- list()
  node.shape <- c(rep("circle", (length(cont.all)) ), rep("box",(length(disc.all)) ))
  names(node.shape) <- c(cont.all, disc.all)
  nAttrs$shape <- node.shape

  if(!is.null(fontsize)){
    nAttrs$fontsize <- rep(fontsize, length(node.names))
    names(nAttrs$fontsize) <- node.names
  }

  max.kl <- max(klds)
  min.kl <- -max(klds)

  if (max.kl==0) {
    fill.post <- rep("white", length(klds.0))
  } else {
    pseudo <- c(klds.0, -klds.0)/max.kl
    # cl.kl <- (klds-min.kl)/(max.kl-min.kl)
    pseudo <- sign(pseudo)*abs(pseudo)^0.5

    rbPal <- colorRampPalette(c('dodgerblue','white','red'))
    fill.post <- rbPal(100)[as.numeric(cut(pseudo,breaks = 100))]
    fill.post <- fill.post[1:length(klds.0)]
  }

  names(fill.post) <- all.active
  # fill.evn <- rep("green", length(vars))
  # names(fill.evn) <- vars

  fill.1 <- rep("khaki", length(vars))
  names(fill.1) <- var.1
  fill.2 <- rep("orange", length(vars))
  names(fill.2) <- var.2
  fill.inter <- rep("gray", length(vars))
  names(fill.inter) <- var.inter

  nAttrs$fillcolor <- c(fill.1, fill.2, fill.inter, fill.post)

  if(length(cont.active)>0){
    par(oma=c(0,0,0,4))
    # Rgraphviz::plot(dag, nodeAttrs=nAttrs, main="")
    graph_plot <- Rgraphviz::layoutGraph(dag, nodeAttrs=nAttrs)
    Rgraphviz::renderGraph(graph_plot)
    par(oma=c( 0,0,0,0.5))
    color.bar(colorRampPalette(c("dodgerblue", "white", "red"))(1000), -signif(max.kl,1), sym=TRUE)
  } else {
    par(oma=c(0,0,0,4))
    # Rgraphviz::plot(dag, nodeAttrs=nAttrs, main="")
    graph_plot <- Rgraphviz::layoutGraph(dag, nodeAttrs=nAttrs)
    Rgraphviz::renderGraph(graph_plot)
    par(oma=c( 0,0,0,0.5))
    color.bar(colorRampPalette(c("white", "red"))(1000), min=0, max=signif(max.kl,1), sym=FALSE)
  }
  return(klds.0)
}


####### function for plotting colorbar

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='', sym=TRUE) {
  scale = (length(lut)-1)/(max-min)
  n.half <- length(lut)/2
  cl <- c()
  for (i in 1:(length(lut)-1)) {
    if (sym) {
      j <- n.half+1+sign(i-n.half)*(abs(i-n.half)/n.half)^0.5*n.half
    } else {
      j <- 1+(i/length(lut))^0.5*length(lut)
    }
    cl[i] <- lut[j]
  }
  image.plot(legend.only=TRUE, col=cl, zlim=c(min(ticks, na.rm=T), max(ticks, na.rm=T)),
             add=TRUE, horizontal=FALSE, legend.shrink=0.3, legend.cex=0.7, legend.lab=title)
}



