
ModelCompileData <- function(data, dag, node.class) {
  
  # data=liver$data; dag=liver$dag; node.class=liver$node.class;
  # data=nci.s; dag=tree.g
  # data=df.2; dag=bn.graph; node.class=node.class;
  
  ###########
  
  # dag, data.frame, node.class
  # data <- dat
  nodes <- names(node.class)# dag@nodes
  dag.graph <- igraph.from.graphNEL(dag)
  
  value.list <- list()
  discrete.nodes <- nodes[node.class]
  continuous.nodes <- nodes[!node.class]
  
  df <- data
  
  ###### convert discrete variables into characters
  
  df[discrete.nodes] <- lapply(df[discrete.nodes], as.character)
  
  #################################################
  
  dat.complete.0 <- df[complete.cases(df),]
  
  ######################
  ## discrete part starts
  ######################
  cpt.pots <- list()
  
  if (length(discrete.nodes)>0){
    
    for (i in 1:length(discrete.nodes)) {
      value.list[[i]] <- sort(unique(dat.complete.0[[discrete.nodes[i]]]))
    }
    names(value.list) <- discrete.nodes
    
    # source("gRain_test_functions.R")
    
    
    
    for (i in 1:length(discrete.nodes)) {
      # i <- 2
      this.node <- discrete.nodes[i]
      this.parents <- names(neighbors(dag.graph, this.node, "in"))
      all.nodes <- c(this.node, this.parents)
      n.nodes <- length(all.nodes)
      this.cpt <- expand.grid(value.list[all.nodes])
      this.df <- df[all.nodes]
      
      tab <- as.data.frame(xtabs(~., this.df))
      pot.joint <- list(cpt=tab[1:n.nodes], prob=tab$Freq/sum(tab$Freq))
      
      if(length(this.parents)==0) {
        pot <- pot.joint
      } else {
        pot <- conditional(pot.joint, this.parents)
      }
      
      cpt.pots[[i]] <- pot
    }
    
    names(cpt.pots) <- discrete.nodes
    
  }
  
  ######################
  ## discrete part ends
  ######################
  
  
  ######################
  ## continuous part starts
  ######################
  
  bags <- list()
  
  if (length(continuous.nodes)>0){
    
    bags <- vector("list", length(continuous.nodes))
    
    for (i in 1:length(continuous.nodes)){
      # i <- 1
      
      # this.bag <- list()
      k <- 1
      
      this.node <- continuous.nodes[i]
      this.parents <- names(neighbors(dag.graph, this.node, "in"))
      
      dat.complete <- df[, c(this.node, this.parents), drop=FALSE]
      dat.complete <- dat.complete[complete.cases(dat.complete), , drop=FALSE]
      
      this.classes <- node.class[this.parents]
      
      ######################
      
      discrete.parents <- this.parents[which(this.classes)]
      continuous.parents <- this.parents[which(!this.classes)]
      
      ######################
      
      if(length(discrete.parents)==0 & length(continuous.parents)==0){
        
        this.bag <- new("LPPotential", 
                        head = this.node,
                        tail = continuous.parents
        )
        
        y <- dat.complete[[this.node]]
        this.bag@const <- mean(y)
        this.bag@variance <- var(y)
        
        bags[[i]] <- this.bag
        next
      }
      
      ###### no discrete parent, but with continuous parents, there is only one linear model
      
      if(length(discrete.parents)==0){
        
        this.bag <- new("LPPotential", 
                        head = this.node,
                        tail = continuous.parents,
                        beta = matrix(NA, nrow=1, ncol=length(continuous.parents))
        )
        
        colnames(this.bag@beta) <- continuous.parents
        
        df.2 <- dat.complete[, c(this.node, continuous.parents), drop=FALSE]
        df.sub <- df.2
        form.str <- paste0(this.node, "~.")
        form <- as.formula(form.str)
        lm.fit <- lm(form, df.sub)
        coefs <- coef(lm.fit)
        
        this.bag@beta[1,] <- coefs[2:length(coefs)]
        this.bag@const <- coefs[1]
        this.bag@variance <- summary(lm.fit)$sigma^2
        
        bags[[i]] <- this.bag
        next
      }
      
      ####### conditions where there are discrete parents
      
      this.disc.vals <- value.list[discrete.parents]
      this.all.combs <- expand.grid(this.disc.vals, stringsAsFactors=FALSE)
      comb.val.list <- apply(this.all.combs, 1, paste0, collapse="%")
      df.1 <- dat.complete[discrete.parents]  ## mark
      df.comb <- apply(df.1, 1, paste0, collapse="%")
      
      ### condition 3, only discrete parents
      
      if(length(continuous.parents)==0){
        this.bag <- new("LPPotential", 
                        head = this.node,
                        config = matrix(NA, nrow=length(comb.val.list), ncol=length(discrete.parents)),
                        beta = matrix(0, nrow=length(comb.val.list), ncol=0)
        )
        
        colnames(this.bag@config) <- discrete.parents
        
        for (j in 1:length(comb.val.list)) {
          # j <- 1
          sub.ind <- which(df.comb==comb.val.list[j])
          df.2 <- dat.complete[[this.node]]
          y <- df.2[sub.ind]
          
          this.bag@config[j,] <- as.vector(this.all.combs[j,], mode="character")
          this.bag@const[j] <- mean(y)
          this.bag@variance[j] <- var(y)
        }
        bags[[i]] <- this.bag
        next
      }
      
      ### condition 4, both discrete and continuous parents
      
      this.bag <- new("LPPotential", 
                      head = this.node,
                      tail = continuous.parents,
                      config = matrix(NA, nrow=length(comb.val.list), ncol=length(discrete.parents)),
                      beta = matrix(NA, nrow=length(comb.val.list), ncol=length(continuous.parents))
                      )
      
      colnames(this.bag@config) <- discrete.parents
      colnames(this.bag@beta) <- continuous.parents
      
      for (j in 1:length(comb.val.list)) {
        # j <- 1
        sub.ind <- which(df.comb==comb.val.list[j])
        df.2 <- dat.complete[ , c(this.node, continuous.parents), drop=FALSE] #######
        df.sub <- df.2[sub.ind,]
        form.str <- paste0(this.node, "~.")
        form <- as.formula(form.str)
        lm.fit <- lm(form, df.sub)
        coefs <- coef(lm.fit)
        
        this.bag@config[j,] <- as.vector(this.all.combs[j,], mode="character")
        this.bag@beta[j,] <- coefs[2:length(coefs)]
        this.bag@const[j] <- coefs[1]
        this.bag@variance[j] <- summary(lm.fit)$sigma^2
      }
      
      bags[[i]] <- this.bag
    }
    
    names(bags) <- continuous.nodes
    
  }
  
  ######################
  ## continuous part ends
  ######################
  
  result <- list(pots = cpt.pots,
                 bags = bags)
  return(result)
  
}

