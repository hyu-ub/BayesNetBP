
###########################################
## Add Evidence to LPPotential
###########################################

substituteEvidence <- function(lppot, vars, vals)
{
  if(!is.subset(vars, lppot@tail)) {
    stop("The variables are not a subset of the lppotential's tail.")
  }
  
  varinds <- c()
  for(i in 1:length(vars)) {
    varinds[i] <- which(lppot@tail==vars[i])
  }
  
  if (length(varinds)==0) {
    return(lppot)
  }
  
  lppot@const <- c(lppot@const + lppot@beta[,vars, drop=FALSE] %*% vals)
  lppot@tail <- lppot@tail[-varinds]
  lppot@beta <- lppot@beta[, -varinds, drop=FALSE]
  return(lppot)
}

###########################################
## Absorb Continuous Evidence
###########################################

PushOperation <- function(tree.push, var, val){
  
  # Step 1
  e.seq <- tree.push@node # the node is arranged in elimination order
  k <- which(e.seq == var)
  if (k!=1){
    for (i in 1:(k-1)) {
      # i <- 2
      this.cluster <- e.seq[i]
      
      if( length(tree.push@lppotential[[this.cluster]])>0 ) {
        this.tail <- tree.push@lppotential[[this.cluster]][[1]]@tail
      } else {
        this.tail <- NA
      }
      
      if (var %in% this.tail) {
        this.pot <- tree.push@lppotential[[this.cluster]][[1]]
        this.pot <- substituteEvidence(this.pot, var, val)
        tree.push@lppotential[[this.cluster]][[1]] <- this.pot
      }
    }
  }
  
  # Step 2
  tree.push@postbag[[var]] <- tree.push@lppotential[[var]]
  tree.push@lppotential[[var]] <- list()
  tree.push@activeflag[[var]] <- FALSE
  
  # Step 3
  this.par <- tree.push@parent[[var]]
  this.var <- var
  
  while(!is.na(this.par) && !tree.push@cluster.class[[this.par]]){
    
    tree.push@postbag[[this.par]] <- tree.push@postbag[[this.var]]
    
    ## Check if it is necessary to perform exchange operation
    flag <- tree.push@activeflag[[this.par]]
    if (length(tree.push@lppotential[[this.par]])==0){
      flag <- FALSE
    } else {
      ## this check might be redundant, as it is also checked in Exchange function
      lp.head <- tree.push@lppotential[[this.par]][[1]]@head
      postbag.tail <- tree.push@postbag[[this.par]][[1]]@tail
      if( !lp.head %in% postbag.tail){
        flag <- FALSE
      }
    }
    ##
    
    if (flag) {
      newBag <- Exchange(tree.push@postbag[[this.par]][[1]], tree.push@lppotential[[this.par]][[1]])
      tree.push@postbag[[this.par]][[1]] <- newBag$postbag
      tree.push@lppotential[[this.par]][[1]] <- newBag$lppotential
      tree.push@lppotential[[this.par]][[1]] <- substituteEvidence(tree.push@lppotential[[this.par]][[1]], var, val)
    }
    
    tree.push@postbag[[this.var]] <- list()
    this.var <- this.par
    this.par <- tree.push@parent[[this.var]]
  }
  
  # Step 4
  
  if (is.na(this.par)) {
    return(tree.push)
  } else {
    this.pot <- tree.push@postbag[[this.var]][[1]]
    likelihood <- dnorm(rep(val, length(this.pot@const)), mean=this.pot@const, sd=sqrt(this.pot@variance))
    pot <- list(cpt=this.pot@config, prob=likelihood)
    tree.push@cpt[[this.par]] <- factor.product(pot, tree.push@cpt[[this.par]])
    return(tree.push)
  }
  
}

###########################################
## Absorb Discrete Evidence
###########################################

# tree <- tree.init; var <- "HDL"; val <- "High";

DiscreteEvidence <- function(tree, var, val) {
  
  ## CPT
  
  if (length(tree@cpt) > 0) {
    for (i in 1:length(tree@cpt)) {
      tab <- tree@cpt[[i]]$cpt
      prob <- tree@cpt[[i]]$prob
      
      if(var %in% colnames(tab)){
        k <- which(colnames(tab)==var)
        keep <- which(tab[,k] == val)
        
        tree@cpt[[i]]$cpt <- tab[keep, -k, drop=FALSE]
        tree@cpt[[i]]$prob <- tree@cpt[[i]]$prob[keep]
        tree@cpt[[i]]$prob <- tree@cpt[[i]]$prob/sum(tree@cpt[[i]]$prob)
      }
    }
  }

  ## JPT
  
  if (length(tree@jpt) > 0) {
    for (i in 1:length(tree@jpt)) {
      tab <- tree@jpt[[i]]$cpt
      prob <- tree@jpt[[i]]$prob
      
      if(var %in% colnames(tab)){
        k <- which(colnames(tab)==var)
        keep <- which(tab[,k] == val)
        
        tree@jpt[[i]]$cpt <- tab[keep, -k, drop=FALSE]
        tree@jpt[[i]]$prob <- tree@jpt[[i]]$prob[keep]
        tree@jpt[[i]]$prob <- tree@jpt[[i]]$prob/sum(tree@jpt[[i]]$prob)
      }
    }
  }
  
  ## LPPotential
  
  if (length(tree@lppotential) > 0) {
    for (i in 1:length(tree@lppotential)) {
      tab <- tree@lppotential[[i]][[1]]@config
      
      if(var %in% colnames(tab)){
        k <- which(colnames(tab)==var)
        keep <- which(tab[,k] == val)
        
        tree@lppotential[[i]][[1]]@config <- tab[keep, -k, drop=FALSE]
        tree@lppotential[[i]][[1]]@beta <- tree@lppotential[[i]][[1]]@beta[keep, , drop=FALSE]
        tree@lppotential[[i]][[1]]@const <- tree@lppotential[[i]][[1]]@const[keep]
        tree@lppotential[[i]][[1]]@variance <- tree@lppotential[[i]][[1]]@variance[keep]
      }
    }
  }
    
  ## Postbag
  
  if (length(tree@postbag) > 0) {
    for (i in 1:length(tree@postbag)) {
      if (length(tree@postbag[[i]]) > 0) {
        for (j in 1:length(tree@postbag[[i]])) {
          tab <- tree@postbag[[i]][[j]]@config
          
          if(var %in% colnames(tab)){
            k <- which(colnames(tab)==var)
            keep <- which(tab[,k] == val)
            
            tree@postbag[[i]][[j]]@config <- tab[keep, -k, drop=FALSE]
            tree@postbag[[i]][[j]]@beta <- tree@postbag[[i]][[j]]@beta[keep, , drop=FALSE]
            tree@postbag[[i]][[j]]@const <- tree@postbag[[i]][[j]]@const[keep]
            tree@postbag[[i]][[j]]@variance <- tree@postbag[[i]][[j]]@variance[keep]
          }
        }
      }
    }
  }
  return(tree)
}


###########################################
## Virtual Evidence
###########################################

# tree <- tree.init; var <- vars[1]; val <- values[[1]];

VirtualEvidence <- function(tree, var, val) {
  
  df.temp <- data.frame(names(val))
  names(df.temp) <- var
  lk.pot <- list(cpt = df.temp, prob = val/sum(val))
  
  if (length(tree@cpt) > 0) {
    for (i in 1:length(tree@cpt)) {
      this.pot <- tree@cpt[[i]]
      tab <- this.pot$cpt
      
      if(var %in% colnames(tab)){
        result <- factor.product(this.pot, lk.pot)
        tree@cpt[[i]] <- result
        return(tree)
      }
    }
  } else {
    stop("No discrete variable found.")
    return(tree)
  }
}
