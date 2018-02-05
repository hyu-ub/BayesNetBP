#' Compute signed and symmetric Kullback-Leibler divergence 
#'
#' Compute signed and symmetric Kullback-Leibler divergence of variables over a spectrum of evidence
#'
#' @details Compute signed and symmetric Kullback-Leibler divergence of variables over a spectrum of evidence.
#' The signed and symmetric Kullback-Leibler divergence is also known as Jeffery's signed information (JSI) for
#' continuous variables.
#'
#' @param tree a \code{\linkS4class{ClusterTree}} object
#' @param var0 the variable to have evidence absrobed
#' @param vars the variables to have divergence computed
#' @param seq a \code{vector} of numeric values as the evidences
#' @param pbar \code{logical(1)} whether to show progress bar
#' @param method method for divergence computation:
#' \code{gaussian} for Gaussian approximation, \code{} for Monte Carlo integration
#' @return a \code{data.frame} of the divergence
#' 
#' @author Han Yu
#' 
#' @examples 
#' \dontrun{
#' data(liver)
#' tree.init.p <- Initializer(dag=liver$dag, data=liver$data, 
#'                            node.class=liver$node.class, 
#'                            propagate = TRUE)
#' klds <- ComputeKLDs(tree=tree.init.p, var0="Nr1i3", 
#'                     vars=setdiff(tree.init.p@node, "Nr1i3"),
#'                     seq=seq(-3,3,0.5))
#' head(klds)
#' }
#' @export

ComputeKLDs <- function(tree, var0, vars, seq, pbar=TRUE, method = "gaussian") {
  
  # cat(method, "\n")
  
  node.class <- tree@node.class
  tree.graph <- tree@graph$tree
  
  x.seq <- seq
  posteriors.1 <- Marginals(tree, vars)
  n.v <- length(vars)
  klds <- matrix(NA, nrow=length(x.seq), ncol=n.v)
  sys <- Sys.info()[1]
  
  if(pbar){
    if(sys=="Windows") {
      pb <- winProgressBar(title = "Computing divergence", min = 0,
                           max = length(x.seq), width = 300)
    } else {
      pb <- txtProgressBar(title = "Computing divergence", 
                           min = 0, max = length(x.seq), style = 3)
    }
  }
  
  for (i in 1:length(x.seq)) {
    object.2 <- AbsorbEvidence(tree, var0, list(x.seq[i]))
    posteriors.2 <- Marginals(object.2, vars)
    for(j in 1:n.v){
      klds[i,j] <- SymmetricKLD(posteriors.1$marginals[[j]], 
                                posteriors.2$marginals[[j]], 
                                discrete = node.class[vars[j]],
                                method = method) ######
    }
    
    if(pbar){
      if(sys=="Windows") {
        setWinProgressBar(pb, i, title=paste("Computing divergence:", round(i/length(x.seq)*100, 0), "% complete"))
      } else {
        setTxtProgressBar(pb, i, title=paste("Computing divergence:", round(i/length(x.seq)*100, 0), "% complete"))
      }
    }
    
  }
  if(pbar){close(pb)}
  
  df <- data.frame(x.seq, klds)
  colnames(df) <- c("x", vars)
  return(df)
}