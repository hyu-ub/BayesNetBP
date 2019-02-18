
setClass("LPPotential",
         slots = list(head = "character", 
                      tail = "character",
                      config = "matrix",
                      beta = "matrix", # each row for a configuration
                                       # each column for a tail variable
                      const = "numeric",
                      variance = "numeric"
         )
)
