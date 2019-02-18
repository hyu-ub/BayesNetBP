
PlotPosteriorDiscrete <- function(posteriors, group=NULL) {
  
  if(is.null(group)){
    nms <- names(posteriors)
  } else {
    nms <- group
  }
  
  ncl <- length(posteriors[[1]])
  if (ncl<3) ncl <- 3
  
  barplot(t(t(posteriors[[1]])), width=0.5, space=0.2, main=nms,
          col=brewer.pal(ncl, "Blues"), xlab=names(posteriors),
          beside=TRUE, names.arg=names(posteriors[[1]]))
  
}
