
setGeneric("getValue",
           def=function(nodeObj){
             standardGeneric("getValue")
           })

# #' Get the value of a node
# #' @param nodeObj The Node object
# #' @return the value of the Node
# #' @export
setMethod(f = "getValue",
          signature = "Node",
          definition = function(nodeObj){
            return(nodeObj@value)
          })
