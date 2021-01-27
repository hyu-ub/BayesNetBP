
##################################################
## Obtain location of a locus
##################################################

loci.loc <- function(locus) {
  split <- strsplit(locus, split='_')
  chr <- c()
  location <- c()
  for (i in 1:length(split)) {
    chr[i] <- substring(split[[i]][1], 4)
    location[i] <- split[[i]][2]
  }

  result <- data.frame(chr, location)
  result$chr <- as.character(result$chr)
  result$location <- as.numeric(as.character(result$location))

  return(result)
}


#####################################
## function for assigning universe
#####################################

assignUniverse <- function(dag, universes, nodes){
  # universes <- cs.2
  assignment <- list()
  node.names <- dag@nodes
  assigned <- c()
  assigned <- setdiff(dag@nodes, nodes)
  dag.graph <- igraph.from.graphNEL(dag)

  i <- 1
  for (universe in universes){
    temp <- c()
    for (node in universe){
      if (node %in% assigned) next
      node.parents <- names(neighbors(dag.graph, node, mode="in"))
      if (length(node.parents)==0) {
        temp <- c(temp, node)
        assigned <- c(assigned, node)
        next
      }
      if (prod(node.parents %in% universe) == 1) {
        temp <- c(temp, node)
        assigned <- c(assigned, node)
      }
    }
    assignment[[i]] <- temp
    i <- i+1
  }
  names(assignment) <- names(universes)
  return(assignment)
}


############################################################
## a helper function for checking if x is a subset of y ####
############################################################
is.subset <- function(x,y){
  if (length(setdiff(x,y))==0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

############################################################
## Extract info from qtl fit results
############################################################
#' @importFrom graph nodes
#' @importFrom qtlnet loci.qtlnet
#' @importFrom qtl scanone
#' @importFrom igraph is_dag
extractQTL <- function(qtl.fit) {

  qtl.graph <- igraph_from_qtlnet(qtl.fit) # igraph.qtlnet(qtl.fit)
  qtl <- qtl.fit$cross

  dag <- igraph.to.graphNEL(qtl.graph)

  if (!is_dag(qtl.graph)) stop("Graph is not a DAG.")

  graph::nodes(dag) <- gsub("@", "_", graph::nodes(dag))
  node.names <- graph::nodes(dag)

  # pheno <- qtl$pheno[,pheno.cols]
  pheno <- qtl$pheno

  loci <- qtlnet::loci.qtlnet(qtl.fit)
  locus <- unique(unlist(loci))

  locus <- gsub("@", "_", locus)
  qtl.df <- loci.loc(locus)

  discrete.nodes <- locus
  continuous.nodes <- names(pheno)

  node.class <- node.names %in% discrete.nodes
  names(node.class) <- node.names

  geno.list <- list()
  markers <- c()

  for(i in 1:nrow(qtl.df)) {
    markers[i] <- find.marker(qtl, qtl.df$chr[i], qtl.df$location[i])
    geno.list[[i]] <-
      data.frame(qtl$geno[[qtl.df$chr[i]]]$data)[[markers[i]]]
  }

  geno <- matrix(unlist(geno.list), byrow=FALSE, ncol=length(geno.list))

  colnames(geno) <- locus

  dat <- data.frame(geno, pheno)

  result <- list(data=dat,
                 dag=dag,
                 node.class=node.class)
  return(result)
}

###############

igraph_from_qtlnet <- function(qtl.fit) {
  qtl_sum <- summary(qtl.fit)[[2]][, 1:2]
  loci_list <- qtlnet::loci.qtlnet(qtl.fit)
  df_loci <- data.frame(cause = unlist(loci_list), effect = names(loci_list))
  rownames(df_loci) <- NULL
  df_edge_list <- rbind(df_loci, qtl_sum)
  igraph_output <- igraph::graph.data.frame(df_edge_list)
  return(igraph_output)
}

