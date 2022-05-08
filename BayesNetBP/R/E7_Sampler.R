#' Sampling from the Bayesian network
#'
#' Sampling from the joint distribution of all applicable nodes in the Bayesian network.
#'
#' @param tree a \code{\linkS4class{ClusterTree}} object
#' @param n a \code{integer} number of observations to generate
#' @return a \code{dataframe} of generated data
#'
#' @author Han Yu
#'
#' @references Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks.
#' Journal of Machine Learning Research, 6(Sep), 1517-1550. \cr
#' \cr
#' Yu H, Moharil J, Blair RH (2020). BayesNetBP: An R Package for Probabilistic Reasoning in Bayesian
#' Networks. Journal of Statistical Software, 94(3), 1-31. <doi:10.18637/jss.v094.i03>.
#'
#' @import doBy
#' @importFrom graph nodes
#' @importFrom igraph neighbors
#' @importFrom methods new
#' @examples
#'
#' data(toytree)
#' Sampler(tree = toytree, n = 10)
#'
#' @export

Sampler <- function(tree, n) {
  # tree <- tree.post.2; n <- 100
  abd <- tree@absorbed.variables
  discrete.nodes <- names(tree@node.class)[tree@node.class]
  continuous.nodes <- names(tree@node.class)[!tree@node.class]
  disc.v <- setdiff(discrete.nodes, abd)
  cont.v <- setdiff(continuous.nodes, abd)

  ## Special case, all discrete nodes observed

  if (length(disc.v) == 0) {
    cont.g <- data.frame()
    for (i in 1:n){
      vec.g <- continuous.single.sampler.special(tree, cont.v)
      cont.g <- rbind(cont.g, vec.g)
    }
    colnames(cont.g) <- cont.v
    rownames(cont.g) <- NULL
    return(cont.g)
  }

  ###########################################

  disc.jd <- FactorQuery(tree, vars = disc.v, mode = "joint")
  cnts <- rmultinom(n = 1, size = n, prob = disc.jd$prob)
  config.tab <- disc.jd[, 1:(ncol(disc.jd)-1)]
  config.tab <- data.frame(lapply(config.tab, as.character), stringsAsFactors=FALSE)

  cont.g <- data.frame()
  for (i in 1:nrow(config.tab)) {
    if (cnts[i] == 0) {
      next
    }
    # i <- 1
    this.config <- unlist(config.tab[i, , drop = TRUE])
    ## generate continuous variables
    for (j in 1:cnts[i]) {
      vec.g <- continuous.single.sampler(tree, cont.v, this.config)
      cont.g <- rbind(cont.g, vec.g)
    }
  }
  colnames(cont.g) <- cont.v
  disc.g <- config.tab[rep(1:nrow(config.tab), cnts), ]
  colnames(disc.g) <- colnames(config.tab)
  generated <- cbind(disc.g, cont.g)
  rownames(generated) <- NULL
  return(generated)
}

######

compatible <- function(config.1, config.2) {

  var.1 <- names(config.1)
  var.2 <- names(config.2)

  var.b <- intersect(var.1, var.2)

  # no_intersection <- TRUE
  # tryCatch(
  #   {
  #     capture.output(graphNEL(nodes = c(var.1, var.2))) # throws error if duplicate nodes
  #   },
  #   error = function(e){
  #     no_intersection <- FALSE
  #   }, warning = function(e){
  #     no_intersection <- FALSE
  #   })
  #
  # if(no_intersection){
  #   return(TRUE)
  # }


  if (length(var.b)==0) {
    return(TRUE)
  }

  config.sub.1 <- config.1[var.b]
  config.sub.2 <- config.2[var.b]
  return(identical(config.sub.1, config.sub.2))
}

######################################

# new version
continuous.single.sampler <- function(tree, cont.v, this.config) {
  x.cont <- rep(NA, length(cont.v))
  names(x.cont) <- cont.v
  x.gen <- c()
  ## "protime", "ast", "alk", "trig", "copper", "chol", "albumin", "bili"
  for (nd in rev(cont.v)) {
    this.pot <- tree@lppotential[[nd]][[1]]

    if(ncol(this.pot@config) == 0) {
      selectedConfig <- 1
    } else {
      same_named_values <- this.config[intersect(colnames(this.pot@config), names(this.config))] # new implementation
      selectedConfig <- which(apply(this.pot@config, 1, function(x) identical(x, same_named_values)))

      # print(c("CONFIG:", selectedConfig))

      # compat <- apply(this.pot@config, 1, compatible, config.2 = this.config)
      # selectedConfig <- which(compat)
    }

    if(length(selectedConfig) > 1){
      warning("More than one configuration selected!")
    }

    if(ncol(this.pot@beta) == 0) {
      mu <- this.pot@const[selectedConfig]
    } else {
      this.beta <- this.pot@beta
      beta.var <- colnames(this.beta)
      var.g <- intersect(x.gen, beta.var)
      betas <- this.beta[selectedConfig, var.g]
      mu <- this.pot@const[selectedConfig] + sum(betas * x.cont[var.g])
    }

    sd <- sqrt(this.pot@variance[selectedConfig])

    x.cont[nd] <- rnorm(1, mean = mu, sd = sd)
    x.gen <- c(x.gen, nd)
  }
  return(x.cont)
}

######################################
## Sampler for no discrete variables
######################################

continuous.single.sampler.special <- function(tree, cont.v) {

  x.cont <- rep(NA, length(cont.v))
  names(x.cont) <- cont.v
  x.gen <- c()

  for (nd in rev(cont.v)) {
    this.pot <- tree@lppotential[[nd]][[1]]

    if(ncol(this.pot@beta) == 0) {
      mu <- this.pot@const[1]
    } else {
      this.beta <- this.pot@beta
      beta.var <- colnames(this.beta)
      var.g <- intersect(x.gen, beta.var)
      betas <- this.beta[1, var.g]
      mu <- this.pot@const[1] + sum(betas * x.cont[var.g])
    }

    sd <- sqrt(this.pot@variance[1])

    x.cont[nd] <- rnorm(1, mean = mu, sd = sd)
    x.gen <- c(x.gen, nd)
  }

  return(x.cont)
}

