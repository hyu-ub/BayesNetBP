
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

###########################################
## Divergence calculation
###########################################

SymmetricKLD <- function(post1, post2, discrete) {
  
  if (!discrete) {
    return(SymKLD.continuous(post1, post2))
  } else {
    return(SymKLD.discrete(post1, post2))
  }
  
}

###########################################
## Divergence for continuous node
###########################################

density.generator <- function(g,x) {
  return(g[1]*dnorm(x, mean=g[2], sd=sqrt(g[3]) ))
}

SymKLD.continuous <- function(post1, post2) {
  step <- 0.01
  x <- seq(-20,20,by=step)
  
  f1 <- apply(post1, 1, density.generator, x)
  y1 <- rowSums(f1)
  
  f2 <- apply(post2, 1, density.generator, x)
  y2 <- rowSums(f2)
  
  m1 <- sum(step*x*y1)
  m2 <- sum(step*x*y2)
  
  kld1 <- sum(step*y1*log(y1/y2))
  kld2 <- sum(step*y2*log(y2/y1))
  return(0.5*(kld1+kld2)*sign(m2-m1))
}

###########################################
## Divergence for discrete node
###########################################

SymKLD.discrete <- function(p1, p2) {
  kld1 <- sum(p1*log(p1/p2))
  kld2 <- sum(p2*log(p2/p1))
  return(0.5*(kld1+kld2))
}