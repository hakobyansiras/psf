Node <- setClass(
  "Node",
  slots = c(
    ID = "numeric",
    value = "numeric",
    name = "character",
    level = "numeric",
    signals = "list"),

  prototype = list(
    ID = -1,
    value = 1,
    name = "",
    level = 0,
    signals = list()),

  validity = function(object)
  {
    if(object@level < 0)
      return("Negative Node level")
    if(object@ID < 0)
      return("Negative Node ID")
    return(T)
  }

)
