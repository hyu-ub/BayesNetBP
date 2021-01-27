
PlotPosteriorContinuous <- function(posteriors, main="", groups=NULL, col=NULL) {
  n <- length(posteriors)
  xys <- PlotPosteriorXY(posteriors)

  if(is.null(groups)) {
    groups <- names(posteriors)
  }

  if(length(groups)!=n){
    groups <- names(posteriors)
    warning("Number of groups does not match number of marginals.")
  }

  x <- c()
  y <- c()
  for(i in 1:n){
    x <- c(x, xys[[i]][[1]])
    y <- c(y, xys[[i]][[2]])
  }

  x.up <- max(x)
  x.low <- min(x)
  y.up <- max(y)
  y.low <- min(y)

  if(is.null(col)) {
    colors <- 1:n
  } else {
    colors <- col
  }

  # series <- names(posteriors)

  graphics::plot.default(xys[[1]][[1]], xys[[1]][[2]],
       xlim=c(x.low, x.up), ylim=c(y.low, y.up), type="l",
       main=main, xlab="", ylab="Density", col=colors)

  if(n>1) {
    for(i in 2:n) {
      lines(xys[[i]][[1]], xys[[i]][[2]], col=colors[i])
    }
  }

  legend("topright", legend=groups, fill=colors, bty="n")
}

###########################################
## Genereate data for plotting
###########################################

PlotPosteriorXY <- function(posteriors) {
  n <- length(posteriors)
  result <- vector("list", n)
  for (j in 1:n){
    msd <- lapply(posteriors, MeanSD)[[j]]
    lower <- msd[1]-5*msd[2]
    upper <- msd[1]+5*msd[2]
    x <- seq(lower, upper, by=0.001)
    y <- rep(0, length(x))
    post.df <- posteriors[[j]]
    for (i in 1:nrow(post.df)){
      y <- y + post.df[i,1]*dnorm(x,
                                  mean=post.df[i,2],
                                  sd=sqrt(post.df[i,3]))
    }
    this.xy <- list(x,y)
    result[[j]] <- this.xy
  }
  names(result) <- names(posteriors)
  return(result)
}

###########################################
## Means and SDs of Mixture Gaussians
###########################################

MeanSD <- function(df) {
  Mean <- sum(df$prob*df$mu)
  SD <- sqrt( sum(df$prob*df$sd^2) +
                sum(df$prob*df$mu^2) -
                (sum(df$prob*df$mu))^2 )
  return(c(Mean=Mean, SD=SD))
}






