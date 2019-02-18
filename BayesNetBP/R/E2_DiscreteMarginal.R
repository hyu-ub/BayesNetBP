############################################################
## Function for getting marginal of a discrete node
############################################################

DiscreteMarginal <- function(tree, var){
  if (length(tree@jpt)==0) stop("Joint distribution table not found. 
                                No discrete nodes in model or the ClusterTree is not propagated.")
  for (i in 1:length(tree@jpt)) {
    j.pot <- tree@jpt[[i]]
    tab <- j.pot$cpt
    if (var %in% colnames(tab)) {
      pot <- marginalize.discrete(j.pot, var)
      v.name <- as.vector(t(pot$cpt))
      result <- pot$prob
      names(result) <- v.name
      return(result)
    }
  }
}