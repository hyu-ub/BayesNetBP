library("shiny")
# Define UI for application that draws a histogram

shinyUI(
  fluidPage(h4("BayesNetBP"),
            fluidRow(
              column(4,
                      h4("Load Model"),
                      # actionButton("compile", label = "Compile"),
                      fileInput("file1", "Choose your file"),
                      actionButton("update", label = "Update"),
                      
                     hr(),
                     h4("Subgraph"),
                     actionButton("sg_exp", label = "Expand"),
                     actionButton("sg_add", label = "Add to list"),
                     actionButton("sg_sub", label = "Subset"),
                     actionButton("sg_reset", label = "Reset"),
                     radioButtons("direct","Expand direction", 
                                  c("Up" = "in", "Down" = "out", "Both"="all"), 
                                  inline=TRUE),
                     #verbatimTextOutput("added"),
                     
                     selectInput("layout", 
                                 label = "Choose graph layout",
                                 choices = c("dagre", "cose")),
                     
                     helpText("Selected Node"),
                     verbatimTextOutput("clickedNode")
                     
                      
              ),
              
              column(4,
                     h4("Fixed Evidence"),
                     selectInput("var", 
                                 label = "Choose variables to set evidence",
                                 choices = tree.init.p@node ),
                     
                     selectInput("getvalue", 
                                 label = "Select value for discrete variable",
                                 choices = c() ),
                     
                     textInput("evidence", "Enter value for continuous variable", 
                               value = "", width = NULL, placeholder = NULL),
                     
                     actionButton("add", label = "Add evidence"),
                     actionButton("clear", label = "Clear"),
                     helpText("Evidence to be absorbed:"),
                     verbatimTextOutput("added"),
                     
                     hr(),
                     actionButton("absorb", label = "Absorb"),
                     actionButton("reset", label = "Reset"),
                     actionButton("marg", label="Plot Marginals"),
                     actionButton("post", label = "Shift in marginals")
              
              ),
                     
              column(4,
                     h4("Effects of a spectrum of evidence"),
                     
                     selectInput("abvar", 
                                 label = "Choose the observed variable",
                                 choices = names(tree.init.p@node.class)[!tree.init.p@node.class]),
                     
                     #selectInput("var2", 
                     #            label = "Choose variables to plot",
                     #            choices = tree.init.p$nodes),
                     
                     actionButton("addplot", label = "Add to plot list"),
                     actionButton("addall", label = "Add all"),
                     actionButton("clearpvar", label = "Clear"),
                     actionButton("plotkld", label = "Plot"),
                     
                     helpText("Variables to plot:"),
                     verbatimTextOutput("addplot"),
                     
                     sliderInput("range", "Range:",
                                 min = -10, max = 10, value = c(-3,3)),
                     sliderInput("increment", "Step:",
                                 min = 0, max = 1, value = 0.5, step= 0.1)
                     
              )
            ),     
                      
            fluidRow(
              column(7,
                     rcytoscapejsOutput("plot", height="700px"),
                     plotOutput("plot1", height = 50, width = 50)
              ),
              
              column(5,
                     htmlOutput("plot3")
              )
            )
              
            
)
)

