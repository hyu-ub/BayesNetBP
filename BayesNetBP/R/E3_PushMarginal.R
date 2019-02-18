
############################################################
## Function for getting marginal of a continuous node
############################################################

# tree.push <- tree.post; var <- "Neu1"

PushMarginal <- function(tree.push, var){
  
  ######
  ######
  tree.push@postbag[[var]] <- tree.push@lppotential[[var]]
  
  ######
  ######
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
    }
    
    tree.push@postbag[[this.var]] <- list()
    this.var <- this.par
    this.par <- tree.push@parent[[this.var]]
  }
  
  ######
  ######
  if (is.na(this.par)) {
    # pure continuous scenario
    marg <- data.frame(prob=1, 
                       mu=tree.push@postbag[[this.var]][[1]]@const[1], 
                       sd=sqrt( tree.push@postbag[[this.var]][[1]]@variance[1] ))
    rownames(marg) <- NULL
  } else {
    this.pot <- tree.push@postbag[[this.var]][[1]]
    disc.pars <- colnames(this.pot@config)
    par.marg <- marginalize.discrete(tree.push@jpt[[this.par]], disc.pars)
    ind <- index.generator(this.pot@config, par.marg$cpt)
    ind.1 <- ind[,1]
    ind.2 <- ind[,2]
    
    marg <- data.frame(prob=par.marg$prob[ind.2], 
                       mu=this.pot@const[ind.1], 
                       sd=sqrt( this.pot@variance[ind.1] ))
    rownames(marg) <- NULL
  }
  return(marg)
}
