Edge <- setClass(
  Class = "Edge",

  contains = "Node",

  slots = c(
    edgeType = "character",
    weight = "numeric",
    signal = "numeric",
    rank = "numeric",
    loopCount = "numeric",
    isBackward = "logical",
    sourceNode = "Node",
    targetNode = "Node"),



  prototype = list(
    edgeType = "activation",
    weight = 1,
    signal = 0,
    rank = 0,
    loopCount = 0,
    isBackward = F,
    sourceNode = NULL,
    targetNode = NULL),

  validity = function(object)
  {
    if(is.null(object@sourceNode))
      return("NULL sourceNode")
    if(is.null(object@targetNode))
      return("NULL targetNode")
    return(T)
  }

)
