
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

SymmetricKLD <- function(post1, post2, discrete, method = "gaussian", epsilon = 10^-6) {
  if (!discrete) {
    if (method == "mc"){
      return(SymKLD.continuous.mc(post1, post2))
    } else {
      return(SymKLD.continuous.gaussian(post1, post2))
    }
    
  } else {
    return(SymKLD.discrete(post1, post2, epsilon = epsilon))
  }
}

###########################################
## Divergence for continuous node
###########################################

density.generator <- function(g,x) {
  return(g[1]*dnorm(x, mean=g[2], sd=g[3] ))
}

# calclulate KL divergence by Monte Carlo integration

SymKLD.continuous.mc <- function(post.1, post.2) {
  
  n.sub.1 <- nrow(post.1)
  n.sub.2 <- nrow(post.2)
  
  ######
  
  if (n.sub.1 == 1 && n.sub.2 == 1) {
    mu1 <- post.1$mu[1]
    mu2 <- post.2$mu[1]
    sigma1 <- post.1$sd[1]
    sigma2 <- post.2$sd[1]
    kld.1 <- log(sigma2/sigma1) + (sigma1^2 + (mu1 - mu2)^2)/(2*sigma2^2) - 0.5
    kld.2 <- log(sigma1/sigma2) + (sigma2^2 + (mu2 - mu1)^2)/(2*sigma1^2) - 0.5
    kld <- sign(mu2-mu1) * (kld.1 + kld.2)/2
    return(kld)
  }
  
  ######
  
  cnt <- rmultinom(1, 10000, prob = post.1$prob)
  x <- c()
  for (i in 1:n.sub.1) {
    this.x <- rnorm(cnt[i], post.1$mu[i], post.1$sd[i])
    x <- c(x, this.x)
  }
  
  f1 <- apply(post.1, 1, density.generator, x)
  f2 <- apply(post.2, 1, density.generator, x)
  y1 <- rowSums(f1)
  y2 <- rowSums(f2)
  
  kld.1 <- mean(log(y1/y2))
  
  #########
  
  cnt <- rmultinom(1, 10000, prob = post.2$prob)
  x <- c()
  for (i in 1:n.sub.2) {
    this.x <- rnorm(cnt[i], post.2$mu[i], post.2$sd[i])
    x <- c(x, this.x)
  }
  
  f1 <- apply(post.1, 1, density.generator, x)
  f2 <- apply(post.2, 1, density.generator, x)
  y1 <- rowSums(f1)
  y2 <- rowSums(f2)
  
  kld.2 <- mean(log(y2/y1))
  
  mu1 <- sum(post.1$prob * post.1$mu)
  mu2 <- sum(post.2$prob * post.2$mu)
  kld <- sign(mu2-mu1) * (kld.1 + kld.2)/2
  return(kld)
}

# calclulate KL divergence by Gaussian approximation

SymKLD.continuous.gaussian <- function(post.1, post.2) {
  meansd.1 <- MeanSD(post.1)
  meansd.2 <- MeanSD(post.2)
  names(meansd.1) <- names(meansd.2) <- NULL
  
  mu1 <- meansd.1[1]
  mu2 <- meansd.2[1]
  sigma1 <- meansd.1[2]
  sigma2 <- meansd.2[2]
  kld.1 <- log(sigma2/sigma1) + (sigma1^2 + (mu1 - mu2)^2)/(2*sigma2^2) - 0.5
  kld.2 <- log(sigma1/sigma2) + (sigma2^2 + (mu2 - mu1)^2)/(2*sigma1^2) - 0.5
  kld <- sign(mu2-mu1) * (kld.1 + kld.2)/2
  return(kld)
}

##############

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

SymKLD.discrete <- function(p1, p2, epsilon = 10^-6) {
  
  ind_01 <- which(p1 == 0)
  ind_11 <- which(p1 != 0)
  ind_02 <- which(p2 == 0)
  ind_12 <- which(p2 != 0)
  n1 <- length(ind_01)
  n2 <- length(ind_02)
  
  if(n1 > 0) {
    p1[ind_01] <- epsilon
    compensate <- epsilon*n1
    p1[ind_11] <- p1[ind_11] - compensate*p1[ind_11]
  }
  
  if(n2 > 0) {
    p2[ind_02] <- epsilon
    compensate <- epsilon*n2
    p2[ind_12] <- p2[ind_12] - compensate*p2[ind_12]
  }
  
  kld1 <- sum(p1*log(p1/p2))
  kld2 <- sum(p2*log(p2/p1))
  return(0.5*(kld1+kld2))
}