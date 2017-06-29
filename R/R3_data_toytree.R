#' A synthetic toy dataset
#' 
#' This is a synthetic dataset with five discrete and four continuouos
#' variables. The local models of this dataset were randomly generated.
#' Specifically, the coefficients of linear models were generated from \emph{N(0,1)},
#' the variances of the error terms from \emph{chi-square(1)}, and the conditional 
#' probabilities for discrete factors generated from \emph{uniform(0,1)} followed by 
#' normalization. The data set contains a propagated \code{clustertree} object, 
#' which is ready for evidence absorption and making queries.
#' 
#' @format The data set contains a propagated \code{clustertree} object \code{toytree}, 
#' which is ready for evidence absorption and making queries.
#' 
#' @usage data(toytree)
#' @name toytree
#' @docType data
#' 
NULL