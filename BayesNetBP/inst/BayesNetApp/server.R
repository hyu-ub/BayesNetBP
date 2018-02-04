
shinyServer(function(input, output, session) {
  
  ##########################################
  ## Initialize
  ##########################################
  # tree.init.p <- toytree
  dag <- tree.init.p@graph$dag
  node.class <- tree.init.p@node.class
  
  ## Reactive values
  
  v <- reactiveValues(absorb =FALSE, 
                      reset=FALSE, 
                      post=FALSE, 
                      marg=FALSE,
                      plotkld=FALSE,
                      tree.init=tree.init.p, 
                      tree.post=tree.init.p,
                      vars=c(), 
                      values=list(), 
                      pvars=c(), 
                      df.list=list(), 
                      nodecolor=c(),
                      href=c(),
                      subset=tree.init.p@node,
                      # subset=c("A","B","C","D","E"),
                      pre_sub=c(),
                      sub_sel=c(),
                      cltemp=c(),
                      flag=FALSE
                      )
  
  ##########################################
  ## Initial plot setting
  ##########################################
  
  edgeData <- data.frame(matrix(unlist(gRbase::edgeList(dag)), ncol=2, byrow=TRUE))
  names(edgeData) <- c("source", "target")
  name <- id <- tree.init.p@node
  nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
  nodeData$color <- rep("dodgerblue", nrow(nodeData)) 
  # nodeData$BorderColor <- rep("black", nrow(nodeData))
  nodeData$color[which(node.class)] <- "dodgerblue"
  # nodeData$href <- paste0(nodeData$name)
  
  nodeData$shape <- rep("ellipse", nrow(nodeData))
  nodeData$shape[which(node.class)] <- "roundrectangle"
  
  v$nodecolor <- nodeData$color
  v$href <- nodeData$href
  ##########################################
  ## Reactive control
  ##########################################
  
  ## update file
  
  observeEvent(input$update, {
    #str(load(input$file1$datapath))
    load(input$file1$datapath)
    v$tree.init <- tree.init.p
    v$tree.post <- tree.init.p
    v$vars <- c()
    v$values <- list()
    v$pvars <- c()
    ##
    v$subset=tree.init.p@node
    v$pre_sub=c()
    v$sub_sel=c()
    v$cltemp=c()
    v$href=v$tree.init@node
    
    v$absorb =FALSE 
    v$reset=FALSE
    v$post=FALSE
    v$marg=FALSE
    v$plotkld=FALSE
    
    v$nodecolor <- rep("dodgerblue", length(v$tree.init@node)) 
    
    updateSelectInput(session, inputId="var", label = "Choose variables to set evidence", 
                      choices = v$tree.init@node,
                      selected = v$tree.init@node[1])
    
    updateSelectInput(session, inputId="abvar", label = "Choose the observed variable", 
                      choices = v$tree.init@node[!v$tree.init@node.class],
                      selected = v$tree.init@node[!v$tree.init@node.class][1])
    
    updateSelectInput(session, inputId="var2", label = "Choose variables to plot", 
                      choices = v$tree.init@node,
                      selected = v$tree.init@node[1])
  })
  
  ## Absorb evidence
  observeEvent(input$absorb, {
    v$absorb <- TRUE
    v$reset <- FALSE
    
    text.in <- input$evidence
    text.sep <- strsplit(text.in, split=",")[[1]]
    
    vars <- v$vars
    values <- v$values
    
    if(length(vars)>0){
      v$tree.post <- AbsorbEvidence(v$tree.init, vars, values)
    }
    v$nodecolor[which(v$tree.init@node %in% v$vars)] <- "green"
  })
  
  observeEvent(input$clear, {
    v$vars <- c()
    v$values <- c()
  })
  
  observeEvent(input$marg, {
    v$marg <- TRUE
  })
  
  #####
  observeEvent(input$clearpvar, {
    v$pvars <- c()
    v$plotkld <- FALSE
    v$post <- FALSE
  })
  #####
  
  observeEvent(input$add, {
    
    node.class <- v$tree.init@node.class
    
    if (input$var %in% v$vars) {
      k <- which(v$vars==input$var)
      if (node.class[input$var]) {
        v$values[k] <- input$getvalue
      } else {
        v$values[k] <- as.numeric(input$evidence)
      }
    } else {
      v$vars <- c(v$vars, input$var)
      if (node.class[input$var]) {
        v$values[length(v$values)+1] <- input$getvalue
      } else {
        v$values[length(v$values)+1] <- as.numeric(input$evidence)
      }
    }
  })
  
  observeEvent(input$addplot, {
    v$pvars <- c(v$pvars, input$clickedNode)
  })
  
  #####
  
  observeEvent(input$plotkld, {
    v$plotkld <- TRUE
    df <- ComputeKLDs(tree=v$tree.init, var0=input$abvar, vars=v$pvars, 
                      seq=seq(input$range[1], input$range[2], input$increment), pbar=TRUE)
    v$df.list <- list(df)
  })
  
  ## Reset button
  
  observeEvent(input$reset, {
    v$reset <- TRUE
    v$absorb <- FALSE
    v$post <- FALSE
    v$vars <- c()
    v$values <- c()
    v$tree.post <- v$tree.init
    v$nodecolor <- rep("dodgerblue", length(v$tree.init@node))
    v$href <- v$tree.init@node
  })
  
  ########################################
  ## Subgraph
  ########################################
  
  observeEvent(input$sg_add, {
    
    if(!v$flag) {
      v$cltemp <- v$nodecolor
      v$flag <- TRUE
    }
    
    if(length(v$sub_sel)==0) {
      v$pre_sub <- union(v$pre_sub, input$clickedNode)
    } else {
      v$pre_sub <- union(v$pre_sub, v$sub_sel)
    }
    
    v$sub_sel <- c()
    v$nodecolor[which(v$tree.init@node %in% v$pre_sub)] <- "purple"
  })
  
  observeEvent(input$sg_exp, {
    
    if(!v$flag) {
      v$cltemp <- v$nodecolor
      v$flag <- TRUE
    }
    
    if(length(v$sub_sel)==0){
      v$sub_sel <- input$clickedNode
    }
    
    temp <- c()
    for (i in 1:length(v$sub_sel)){
      this.node <- v$sub_sel[i]
      
      if (input$direct=="in") {
        temp <- union(temp, gRbase::parents(this.node, v$tree.init@graph$dag))
      } else if (input$direct=="out") {
        temp <- union(temp, gRbase::children(this.node, v$tree.init@graph$dag))
      } else {
        temp <- union(temp, gRbase::parents(this.node, v$tree.init@graph$dag))
        temp <- union(temp, gRbase::children(this.node, v$tree.init@graph$dag))
      }
      
      
    }
    
    v$sub_sel <- union(v$sub_sel, temp)
    # cat(v$sub_sel)
    v$nodecolor[which(v$tree.init@node %in% v$sub_sel)] <- "orange"
  })
  
  observeEvent(input$sg_sub, {
    
    if(length(v$pre_sub)==0) return()
    
    v$subset <- v$pre_sub
    v$pre_sub <- c()
    v$sub_sel <- c()
    v$nodecolor <- v$cltemp # rep("dodgerblue", length(v$tree.init$nodes))
    v$flag <- FALSE
  })
  
  observeEvent(input$sg_reset, {
    v$subset <- v$tree.init@node
    v$pre_sub <- c()
    v$sub_sel <- c()
    v$nodecolor <- v$cltemp # rep("dodgerblue", length(v$tree.init$nodes))
    v$flag <- FALSE
  })
  
  ########################################
  ## Plot shift in marginals
  ########################################
  observeEvent(input$post, {
    #v$reset <- TRUE
    v$absorb <- FALSE
    #v$vars <- c()
    #v$values <- c()
    v$post <- TRUE
    
    name <- v$tree.init@node
    
    klds <- PlotCGBN(v$tree.init, v$tree.post, fontsize = 30, plotting =FALSE, pbar=TRUE)
    cls <- color.generator(klds)
    
    klds.ab <- rep("Observed", length(v$vars))
    cls.ab <- rep("green", length(v$vars)) 
    names(cls.ab) <- v$vars
    names(klds.ab) <- v$vars
    
    cls.all <- c(cls, cls.ab)
    cls.all <- cls.all[name]
    
    klds.c <- as.character(signif(klds, 3))
    names(klds.c) <- names(klds)
    
    klds.all <- c(klds.c, klds.ab)
    klds.all <- klds.all[name]
    v$nodecolor <- cls.all
    v$href <- paste0(name, ": JSI: ", klds.all)
    
  })
  
  ## Add all variables to plot
  
  observeEvent(input$addall, {
    var0 <- input$abvar
    v$pvars <- setdiff(v$tree.init@node, var0)
  })
  
  ##########################################
  
  output$added <- renderPrint({
    if (length(v$vars)==0) return ("")
    
    evid <- ""
    for (i in 1:length(v$vars)) {
      this.evid <- paste0(v$vars[i],"=", v$values[i], "; ")
      evid <- paste0(evid, this.evid)
    }
    return(evid)
  })
  
  output$addplot <- renderPrint({
    if (length(v$pvars)==0) return ("")
    return(paste0(v$pvars, collapse=","))
  })
  
  output$clickedNode = renderPrint({
    paste(input$clickedNode)
  })
  
  ########################
  ## Get discrete values
  ########################
  observe({
    x <- input$var
    
    values <- GetValue(v$tree.post, x, message=FALSE)
    
    # Can use character(0) to remove all choices
    if (is.null(values))
      values <- character(0)
    
    # Can also set the label and select items
    updateSelectInput(session, "getvalue",
                      label = "Select value for discrete variable",
                      choices = values,
                      selected = values[1]
    )
  })
  
  
  
  ##########################################
  ## Plot
  ##########################################
  
  output$plot <- renderRcytoscapejs({
    
    #########################################
    
    dag <- v$tree.init@graph$dag
    node.class <- v$tree.init@node.class
    
    edgeData <- data.frame(matrix(unlist(gRbase::edgeList(dag)), ncol=2, byrow=TRUE))
    names(edgeData) <- c("source", "target")
    name <- id <- v$tree.init@node
    nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
    nodeData$color <- rep("dodgerblue", nrow(nodeData)) 
    # nodeData$BorderColor <- rep("black", nrow(nodeData))
    nodeData$color[which(node.class)] <- "dodgerblue"
    # nodeData$href <- paste0(nodeData$name)
    
    nodeData$shape <- rep("ellipse", nrow(nodeData))
    nodeData$shape[which(node.class)] <- "roundrectangle"
    
    nodeData$color <- v$nodecolor
    
    #########################################
    
    # nodeData$color[name %in% v$vars] <- "green"
    nodeData$tooltip <- v$href
    # nodeData$href <- as.character(nodeData$href)
    nodeData$color <- v$nodecolor
    
    nodeData.2 <- subset(nodeData, name %in% v$subset)
    edgeData.2 <- subset(edgeData, (source %in% v$subset) & (target %in% v$subset))
    
    network <- createCytoscapeJsNetwork(nodeData.2, edgeData.2)
    rcytoscapejs(network$nodes, network$edges, showPanzoom=TRUE, 
                 # layout="dagre", 
                 layout=input$layout,
                 highlightConnectedNodes = FALSE)
    # print(nodeData)
  })
  
  output$plot1 <- renderPlot({
    
    if(v$marg){
      v$marg <- FALSE
      
      node.class <- v$tree.init@node.class
      
      x11()
      if(!v$absorb & !v$post){
        marg <- Marginals(v$tree.post, input$clickedNode)
        PlotMarginals(marg)
      } else {
        marg.1 <- Marginals(v$tree.init, input$clickedNode)
        marg.2 <- Marginals(v$tree.post, input$clickedNode)
        
        mgns <- list(marg.1$marginals[[1]], marg.2$marginals[[1]])
        names(mgns) <- c(input$clickedNode, input$clickedNode)
        types <- c(marg.1$types[1], marg.2$types[1])
        marg <- list(marginals=mgns, types=types)
        PlotMarginals(marg, groups=c("Before", "After"))
      }
    }
    
  })
  
  
  output$plot3 <- renderGvis({
    if(v$plotkld){
      Line <- gvisLineChart(v$df.list[[1]],
                            options=list(width=600, height=600))
      return (Line)
    }
  })
  
  
})


color.generator <- function(klds) {
  klds.0 <- klds
  klds <- abs(klds.0)
  max.kl <- max(klds)
  min.kl <- -max(klds)
  
  if (max.kl==0) {
    fill.post <- rep("gray70", length(klds.0))
  } else {
    pseudo <- c(klds.0, -klds.0)/max.kl
    # cl.kl <- (klds-min.kl)/(max.kl-min.kl)
    pseudo <- sign(pseudo)*abs(pseudo)^0.5
    
    rbPal <- colorRampPalette(c('dodgerblue', 'azure2', 'red'))
    fill.post <- rbPal(100)[as.numeric(cut(pseudo,breaks = 100))]
    fill.post <- fill.post[1:length(klds.0)]
  }
  
  names(fill.post) <- names(klds)
  return(fill.post)
  
}


