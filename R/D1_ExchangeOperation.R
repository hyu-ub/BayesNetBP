
###########################################
## Exchange of two LPPotentials
###########################################

Exchange <- function(bag.post, bag.lp) {
  
  # bag.post <- tree.push@postbag[[this.par]][[1]]; bag.lp <- tree.push@lppotential[[this.par]][[1]]
  
  if (!bag.lp@head %in% bag.post@tail) { 
    # if the head of lpp is not in the tail of the postbag potential,
    # there is no need to perform exchange operation
    
    result <- list(postbag = bag.post,
                   lppotential = bag.lp)
    return(result)
    
  }
  
  ## need to deal with special cases where at least one configuation is empty
  
  ind.exp <- index.generator(bag.post@config, bag.lp@config)
  nt <- nrow(ind.exp)
  
  ind.1 <- ind.exp[,1]
  ind.2 <- ind.exp[,2]
  
  beta.1 <- bag.post@beta[ind.1, , drop=FALSE]
  
  ## changed here for pure continuous case
  if (ncol(bag.lp@beta)==0 && nrow(bag.lp@beta)==0) {
    beta.2 <- matrix(0, ncol=0, nrow=1)
  } else {
    beta.2 <- bag.lp@beta[ind.2, , drop=FALSE]
  }
  
  ## processing variables
  
  w <- setdiff(union(bag.post@tail, bag.lp@tail), bag.lp@head) # name of W variables
  b <- beta.1[, bag.lp@head] # vector
  
  rem.1 <- setdiff(w, bag.post@tail)
  rem.2 <- setdiff(w, bag.lp@tail)
  
  if (length(rem.1)==0) {
    a <- beta.1[, w, drop=FALSE]
  } else {
    expand.1 <- matrix(0, nrow=nt, ncol=length(rem.1))
    colnames(expand.1) <- rem.1
    a.0 <- cbind(beta.1, expand.1)
    a <- a.0[, w, drop=FALSE]
  }
  
  if (length(rem.2)==0) {
    if (length(w)==0) {
      c <- beta.2
    } else {
      c <- beta.2[, w, drop=FALSE]
    }
    
  } else {
    expand.2 <- matrix(0, nrow=nt, ncol=length(rem.2))
    colnames(expand.2) <- rem.2
    c.0 <- cbind(beta.2, expand.2)
    c <- c.0[, w, drop=FALSE]
  }
  
  s1 <- bag.post@variance[ind.1]
  s2 <- bag.lp@variance[ind.2]
  const1 <- bag.post@const[ind.1]
  const2 <- bag.lp@const[ind.2]
  
  ## Exchange operation
  
  beta.post <- a + b*c  # matrix
  colnames(beta.post) <- w
  const.post <- const1 + b*const2
  denom <- s.post <- s1 + b^2*s2  # vec
  beta.lp.w <- (c*s1 - a*b*s2)/denom
  beta.z <- b*s2/denom
  beta.lp <- cbind(beta.z, beta.lp.w)
  colnames(beta.lp) <- c(bag.post@head, w)
  const.lp <- (const2*s1 - const1*b*s2)/denom
  s.lp <- s1*s2/denom
  
  ## manipulation of configurations
  
  config.var.1 <- colnames(bag.post@config)
  config.var.2 <- colnames(bag.lp@config)
  config.var.int <- intersect(config.var.1, config.var.2)
  config.rem.1 <- setdiff(config.var.1, config.var.int)
  config.rem.2 <- setdiff(config.var.2, config.var.int)
  
  ###### changed > 1.2.1
  if(ncol(bag.post@config)==0 && ncol(bag.lp@config)==0) {
    config.after <- matrix(0, ncol=0, nrow=0)
  } else if (ncol(bag.post@config)==0) {
    config.after <- bag.lp@config
  } else if (ncol(bag.lp@config)==0) {
    config.after <- bag.post@config
  } else {
    config.after <- cbind(bag.post@config[ind.1, config.rem.1, drop=FALSE], 
                          bag.lp@config[ind.2, config.rem.2, drop=FALSE],
                          bag.lp@config[ind.2, config.var.int, drop=FALSE])
  }
  
  ## wrap up the results
  
  bag.post.after <- new("LPPotential",
                        head = bag.post@head,
                        # tail = colnames(beta.post),
                        config = config.after,
                        beta = beta.post,
                        const = const.post,
                        variance = s.post
  )
  
  ## bag.lp.after@tail will never be empty
  
  if(ncol(beta.post)==0) {
    bag.post.after@tail <- character(0)
  } else {
    bag.post.after@tail <- colnames(beta.post)
  }
  
  bag.lp.after <- new("LPPotential",
                      head = bag.lp@head,
                      # tail = colnames(beta.lp),
                      config = config.after,
                      beta = beta.lp,
                      const = const.lp,
                      variance = s.lp
  )
  
  if(ncol(beta.lp)==0) {
    bag.lp.after@tail <- character(0)
  } else {
    bag.lp.after@tail <- colnames(beta.lp)
  }
  
  result <- list(postbag = bag.post.after,
                 lppotential = bag.lp.after)
  
  return(result)
}