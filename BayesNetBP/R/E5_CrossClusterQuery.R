############################################################
## Function for getting joints for factors across clusters
############################################################

query.ooc <- function(tree, vars){
  
  tree.graph <- igraph.from.graphNEL(tree@graph$tree)
  
  cs.sets <- list()
  
  discrete.clusters <- names(tree@cluster.class)[tree@cluster.class]
  cs.sets <- tree@member[discrete.clusters]
  vars.temp <- vars
  
  cs <- c()
  while(length(vars.temp)>0) {
    maxl <- 0
    
    inter.temp <- c()
    temp <- character(0)
    for (i in 1:length(cs.sets)) {
      
      inter <- intersect(vars.temp, cs.sets[[i]])
      if (length(inter)>=maxl) {
        maxl <- length(inter)
        temp <- discrete.clusters[i]
        inter.temp <- inter
      }
    }
    cs <- c(cs, temp)
    vars.temp <- setdiff(vars.temp, inter.temp)
  }
  
  sub.memb <- c()
  
  for (i in 2:length(cs)) {
    path <- all_simple_paths(tree.graph, cs[1], cs[i], mode="all")[[1]]
    sub.memb <- union(sub.memb, names(path))
  }
  
  sub.graph <- induced_subgraph(tree.graph, sub.memb)
  
  # x11(); plot(sub.graph)
  
  node <- sub.memb[1]
  
  #########################
  #########################
  
  
  ooc <- list(jpt=tree@jpt, sub.graph=sub.graph, active=c(), 
              nom=tree@jpt[[node]], 
              denom=list(cpt=matrix(0,nrow=0,ncol=0), prob=1), 
              cs.sets=cs.sets)
  
  obj <- Distribute.OOC(ooc, node)
  temp.pot <- factor.divide(obj$nom, obj$denom)
  jnt <- marginalize.discrete(temp.pot, vars)
  
  return(jnt)
}

## object.ooc <- ooc

Distribute.OOC <- function(object.ooc, node) {
  ngbs <- neighbors(object.ooc$sub.graph, node, mode = "all")$name
  inactive <- setdiff(ngbs, object.ooc$active)
  object.ooc$active <- c(object.ooc$active, node)
  jpt <- object.ooc$jpt[[node]]
  # cat(cluster@members, "\n")
  if (length(inactive)>0) {
    for (i in 1:length(inactive)) {
      # cat(node, "->", inactive[i], "\n")
      this.jpt <- object.ooc$jpt[[inactive[i]]]
      object.ooc$nom <- factor.product(object.ooc$nom, this.jpt, normalize=FALSE)
      # cat(names(object.ooc$nom$cpt), "\n")
      separator <- intersect(object.ooc$cs.sets[[node]], object.ooc$cs.sets[[inactive[i]]])
      # cat(separator, "\n")
      margin <- marginalize.discrete(jpt, separator)
      object.ooc$denom <- factor.product(object.ooc$denom, margin, normalize=FALSE)
      # object.ooc$joint <- factor.divide(this.cluster@joint, margin)
      object.ooc <- Distribute.OOC(object.ooc, inactive[i])
    }
  }
  return(object.ooc)
}





