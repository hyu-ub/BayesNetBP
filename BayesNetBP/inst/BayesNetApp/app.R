library(shiny)
library(cyjShiny)
library(BayesNetBP)
library(googleVis)
library(igraph)

options(shiny.maxRequestSize = 300*1024^2)

data("toytree")
# `tree` and `subtree` are objects of class ClusterTree
tree <- toytree

## Helper functions ##################

MeanSD <- function(df) {
  Mean <- sum(df$prob*df$mu)
  SD <- sqrt( sum(df$prob*df$sd^2) +
                sum(df$prob*df$mu^2) -
                (sum(df$prob*df$mu))^2 )
  return(c(Mean=Mean, SD=SD))
}

#************************************************************
ui <- fluidPage(
  titlePanel("BayesNetBP"),

  fluidRow(
    column(4,
           h4("Load Model"),
           # actionButton("compile", label = "Compile"),
           fileInput("file1", "Select your file"),
           actionButton("update", label = "Update"),
           actionButton("sg_fit", label = "Fit"),
           actionButton("sg_fitsel", label = "Fit Selected"),
           h4("Subgraph"),
           actionButton("sg_exp", label = "Expand"),
           actionButton("sg_sub", label = "Subset"),
           actionButton("sg_reset", label = "Reset"),
           h4("Set Evidence"),
           selectInput("var",
                       label = "Variable to set evidence",
                       choices = sort(tree@node) )
    ),

    # column(2,),

    column(4,
           h4("Evidence value"),
           selectInput("getvalue",
                       label = "Select value (discrete)",
                       choices = c() ),
           textInput("evidence", "Enter value (continuous)",
                     value = "", width = NULL, placeholder = NULL),
           h4("Observe and Query"),
           actionButton("add", label = "Add evidence"),
           actionButton("clear", label = "Clear"),
           actionButton("absorb", label = "Absorb"),
           helpText("Evidence to absorb:"),
           verbatimTextOutput("added"),
           actionButton("marg", label= "Marginal of Selected"),
           actionButton("post", label = "Shift in marginals"),
           actionButton("reset", label = "Reset")

           #
           # selectInput("abvar",
           #             label = "Choose the observed variable",
           #             choices = tree@node),
           #
           # actionButton("addplot", label = "Add to plot list"),
           # actionButton("addall", label = "Add all"),
           # actionButton("clearpvar", label = "Clear"),
           # actionButton("plotkld", label = "Plot"),
           #
           # helpText("Variables to plot:"),
           # verbatimTextOutput("addplot"),
           #
           # sliderInput("range", "Range:",
           #             min = -10, max = 10, value = c(-3,3)),
           # sliderInput("increment", "Step:",
           #             min = 0, max = 1, value = 0.5, step= 0.1)
    ),

    column(4,
           h4("Effects of a spectrum of evidence"),

           selectInput("abvar",
                       label = "Choose the variable",
                       choices = names(tree@node.class)[!tree@node.class]),

           #selectInput("var2",
           #            label = "Choose variables to plot",
           #            choices = tree.init.p$nodes),

           actionButton("sel_obvar", label = "Select observed"),
           actionButton("addplot", label = "Add to plot"),
           actionButton("addall", label = "Add all"),
           actionButton("clearpvar", label = "Clear"),
           actionButton("plotkld", label = "Plot"),

           helpText("Observed variable:"),
           verbatimTextOutput("obvar"),

           helpText("Variables to plot:"),
           verbatimTextOutput("addplot"),

           #sliderInput("range", "Range:", min = -10, max = 10, value = c(-3,3)),
           #sliderInput("increment", "Step:", min = 0, max = 1, value = 0.5, step= 0.1)
           #numericInput("kld_min", "Min", -10, min = NA, max = NA, step = NA, width = 4),
           #numericInput("kld_max", "Max", 10, min = NA, max = NA, step = NA, width = 4),
           #numericInput("kld_step", "Step", 1, min = NA, max = NA, step = NA, width = 4)

           div(style="display: inline-block; vertical-align:top; width: 90px;",numericInput("kld_min", "Min", -10)),
           div(style="display: inline-block; vertical-align:top; width: 90px;",numericInput("kld_max", "Max", 10)),
           div(style="display: inline-block; vertical-align:top; width: 70px;",numericInput("kld_step", "Step", 1))

  ),

  # graph panels, main graph and subgraph
  fluidRow(
    column(8,
           cyjShinyOutput('cyjShiny')
           # h4("View"),
    ),

    column(4,
           htmlOutput("gtable"),
           htmlOutput("plot3")
    ),

    ## spectrum
    column(4,
           htmlOutput("plot4")
    )
    #########
  ),

  # extra white space at the bottom
  mainPanel(
    h3(textOutput("whitespace"))
  )
)
)

#************************************************************
server <- function(input, output, session) {
  v <- reactiveValues(absorb =FALSE,
                      absorbed = rep(0, length(tree@node)),
                      reset=FALSE,
                      post=FALSE,
                      marg=FALSE,
                      plotkld=FALSE,
                      plotmarg=FALSE,
                      tree.init=tree,
                      tree.post=tree,
                      vars=c(),
                      values=list(),
                      pvars=c(),
                      df.list=list(data.frame(), FALSE),
                      df.kld=list(),
                      newNodes=c()
  )

  tbl.nodes <- data.frame(id=sort(tree@node),
                          type=c(1, 0)[tree@node.class[sort(tree@node)] + 1],
                          kld = 0,
                          absorbed = "n",
                          stringsAsFactors=FALSE)
  elist <- as_edgelist(igraph.from.graphNEL(tree@graph$dag))
  tbl.edges <- data.frame(source=elist[,1],
                          target=elist[,2],
                          interaction=rep("phosphorylates", nrow(elist)),
                          stringsAsFactors=FALSE)

  graph.json <- dataFramesToJSON(tbl.edges, tbl.nodes)
  output$cyjShiny <- renderCyjShiny(
    cyjShiny(graph=graph.json, layoutName="dagre", styleFile="biologicalStyle.js")
  )
  # event observation

  observeEvent(input$update, {
    if(!is.null(input$file1) & class(toytree) == "ClusterTree"){
      tree <- get(load(input$file1$datapath))
      tbl.nodes <- data.frame(id=sort(tree@node),
                              type=c(1, 0)[tree@node.class[sort(tree@node)] + 1],
                              kld = 0,
                              absorbed="n",
                              stringsAsFactors=FALSE)
      elist <- as_edgelist(igraph.from.graphNEL(tree@graph$dag))
      tbl.edges <- data.frame(source=elist[,1],
                              target=elist[,2],
                              interaction=rep("phosphorylates", nrow(elist)),
                              stringsAsFactors=FALSE)
      graph.json <- dataFramesToJSON(tbl.edges, tbl.nodes)
      ###

      output$cyjShiny <- renderCyjShiny(
        cyjShiny(graph=graph.json, layoutName="dagre", styleFile="biologicalStyle.js")
      )

      v$tree.init=tree
      v$tree.post=tree

      updateSelectInput(session, "var",
                        label = "Select variable to set evidence",
                        choices = sort(v$tree.init@node),
                        selected = sort(v$tree.init@node)[1]
      )

      updateSelectInput(session, "abvar",
                        label = "Choose the variable",
                        choices = names(v$tree.init@node.class)[!v$tree.init@node.class],
                        selected = sort(v$tree.init@node)[1]
      )

      # reset all variables
      v$absorbed = rep(0, length(tree@node))
      v$vars=c()
      v$values=list()
      v$df.list=list(data.frame(), FALSE) #
      v$df.kld=list() #
      v$pvars=c()
      v$plotmarg <- FALSE
    }

    # probably remove/change this later
    else{
      print("Invalid file")
    }
  })

  ##############################################################

  observeEvent(input$sg_exp, {
    selectFirstNeighbors(session)
  })

  observeEvent(input$sg_sub, {
    invertSelection(session)
    hideSelection(session)
    invertSelection(session)
  })

  observeEvent(input$sg_reset, {
    showAll(session)
  })

  observeEvent(input$sg_fit, {
    fit(session, padding = 20)
  })

  observeEvent(input$sg_fitsel, {
    fitSelected(session, padding = 20)
  })

  ## Select nodes ################################

  #observeEvent(input$getSelectedNodes, {
  # v$vars <- c(v$vars, getSelectedNodes(session))
  # getSelectedNodes(session)
  #})

  #observeEvent(input$click, ignoreInit=TRUE, {
  #  output$selectedNodesDisplay <- renderText({" "})
  #  getSelectedNodes(session)
  #})

  observeEvent(input$selectedNodes, {
    # getSelectedNodes(session)
    v$newNodes <- input$selectedNodes;
    output$selectedNodesDisplay <- renderText({paste(v$newNodes)})
    # print(v$newNodes)

    this.var <- input$selectedNodes[1]
    this.mar.init.temp <- Marginals(v$tree.init, this.var)
    this.mar.post.temp <- Marginals(v$tree.post, this.var)
    this.mar.init <- this.mar.init.temp[[1]][[1]]
    this.mar.post <- this.mar.post.temp[[1]][[1]]

    if (this.mar.init.temp$types[1]) {
      v$df.list <- list(data.frame(var = names(this.mar.init),
                                   Before = this.mar.init,
                                   After = this.mar.post),
                        TRUE)
    } else {
      msd1 <- MeanSD(this.mar.init)
      msd2 <- MeanSD(this.mar.post)
      x.min <- min(msd1[1] - 3*msd1[2], msd2[1] - 3*msd2[2])
      x.max <- max(msd1[1] + 3*msd1[2], msd2[1] + 3*msd2[2])
      x.vec <- seq(x.min, x.max, length.out = 1000)
      v$df.list <- list(data.frame(x = x.vec,
                                   Before = dnorm(x.vec, msd1[1], msd1[2]),
                                   After = dnorm(x.vec, msd2[1], msd2[2])),
                        FALSE)
    }

  })

  ##############################################################

  output$added <- renderPrint({
    if (length(v$vars)==0) return ("")
    evid <- ""
    for (i in 1:length(v$vars)) {
      this.evid <- paste0(v$vars[i],"=", v$values[i], "; ")
      evid <- paste0(evid, this.evid)
    }
    return(evid)
  })

  #output$clickedNode = renderPrint({
  # paste(v$newNodes)
  #})

  ## Add evidence ############################################

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

  observeEvent(input$clear, {
    v$vars <- c()
    v$values <- list()
  })

  ##############################################

  observe({
    x <- input$var
    values <- GetValue(v$tree.post, x, message=FALSE)

    # Can use character(0) to remove all choices
    if (is.null(values))
      values <- character(0)

    # Can also set the label and select items
    updateSelectInput(session, "getvalue",
                      label = "Select value (discrete)",
                      choices = values,
                      selected = values[1]
    )
  })

  ## Absorb evidence #################
  observeEvent(input$absorb, {
    #v$absorb <- TRUE
    #v$reset <- FALSE

    text.in <- input$evidence
    text.sep <- strsplit(text.in, split=",")[[1]]

    vars <- v$vars
    values <- v$values

    if(length(vars)>0){
      v$tree.post <- AbsorbEvidence(v$tree.init, vars, values)
      # v$absorbed[v$vars] <- 1
      setNodeAttributes(session, "absorbed", vars, rep("y", length(vars)))
    }

  })

  ## Shift in marginals ##########

  observeEvent(input$post, {
    #v$reset <- TRUE
    v$absorb <- FALSE
    #v$vars <- c()
    #v$values <- c()
    v$post <- TRUE

    name <- v$tree.init@node

    klds <- PlotCGBN(v$tree.init, v$tree.post, fontsize = 30, plotting =FALSE, pbar=TRUE)
    # cls <- color.generator(klds)
    klds <- klds/(max(abs(klds)))
    kld.nodes <- names(klds)
    names(klds) <- NULL
    setNodeAttributes(session, "kld", kld.nodes, klds)
  })

  ## Get marginal ###################

  observeEvent(input$marg, {
    v$plotmarg <- TRUE
    getSelectedNodes(session)
    # if(length(v$pre_sub)==0) return()
    #this.var <- input$selectedNodes
    #print(this.var)
    #print(Marginals(v$tree.post, this.var)[[1]][[1]])
    # print(this.mar)
    # v$df.list[[1]] <- data.frame(var = names(this.mar), value = this.mar)
  })

  observeEvent(input$reset, {
    v$vars <- c()
    v$values <- c()
    v$tree.post <- v$tree.init
    setNodeAttributes(session, "kld", v$tree.init@node, rep(0, length(v$tree.init@node)))
  })

  ############

  output$plot3 <- renderGvis({
    if(v$plotmarg){

      if(ncol(v$df.list[[1]]) == 0) return() ##

      if(v$df.list[[2]]) {
        glplot <- gvisBarChart(v$df.list[[1]], options=list(width=400, height=400))
      } else {
        glplot <- gvisLineChart(v$df.list[[1]], options=list(width=400, height=400))
      }
      return (glplot)
    }
  })

  output$gtable <- renderGvis({
    if(v$plotmarg){
      if(v$df.list[[2]]) {
        gtable <- gvisTable(v$df.list[[1]])
      } else {
        gtable <- NULL
      }
      return (gtable)
    }
  })

  #########

  observeEvent(input$sel_obvar, {
    v$obvar <- input$abvar
  })

  output$obvar <- renderPrint({
    if (length(v$obvar)==0) return ("")
    return(v$obvar)
  })

  observeEvent(input$addplot, {
    v$pvars <- c(v$pvars, input$abvar)
  })

  observeEvent(input$clearpvar, {
    v$pvars <- c()
  })

  observeEvent(input$addall, {
    v$pvars <- setdiff(v$tree.init@node, v$obvar)
  })

  observeEvent(input$plotkld, {
    if(length(setdiff(v$pvars, v$obvar)) > 0) {
      v$plotkld <- TRUE
      df <- ComputeKLDs(tree=v$tree.init, var0=v$obvar, vars=setdiff(v$pvars, v$obvar),
                        seq=seq(input$kld_min, input$kld_max, input$kld_step), pbar=TRUE)
      v$df.kld <- list(df)
    }
  })

  output$addplot <- renderPrint({
    if (length(v$pvars)==0) return ("")
    return(paste0(v$pvars, collapse=","))
  })

  output$plot4 <- renderGvis({
    if(v$plotkld){
      Line <- gvisLineChart(v$df.kld[[1]], options=list(width=500, height=500))
      return (Line)
    }
  })

}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
