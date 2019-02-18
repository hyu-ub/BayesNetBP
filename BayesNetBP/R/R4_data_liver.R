#' Mus Musculus HDL QTL data from Leduc et. al. (2012)
#' 
#' Liver QTL data was obtained from a F2 inner-cross between inbred MRL/MpJ and SM/J 
#' strains of mice.
#' 
#' @name liver
#' @usage 
#' data(liver)
#' @format The data set \code{liver} contains three objects: the data, a learned Bayesian network structure 
#' and \code{vector} specifying node type. The fields are described as follows:  
#' 
#' \describe{
#'  \item{\code{data}}{a \code{data.frame} object that contains 280 samples (rows) and 15 variables: genotype data 
#'  (genotype states at 5 SNP markers) and phenotype data (HDL levels and normalized expression values of 10 genes).  
#'  Three of these phenotypes are dichotomized, including Cyp2b10, Spgl1 and HDL.  Genotypes and dichotomized phenotypes 
#'  are of class \code{factor} and continuous phenotypes are of class \code{numeric}.}
#'  \item{\code{dag}}{a \code{graphNEL} object, which is the network structure learned by \code{qtlnet} package.}
#'  \item{\code{node.class}}{a named \code{vector} of \code{logical} values indicating whether each node is discrete.}
#' }
#' 
#' @docType data
#' 
#' @references Leduc MS, Blair RH, Verdugo RA, Tsaih SW, Walsh K, Churchill GA, Paigen B.(2012). 
#' "Using bioinformatics and systems genetics to dissect HDL-cholesterol genetics in an MRL/MpJ 
#' x SM/J intercross." J Lipid Res., 6, 1163-75.

NULL