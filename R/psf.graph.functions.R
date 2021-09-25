out.edges <- function(g, from){
  if(!(from %in% g@nodes)){
    cat("no node ", from, " in the graph\n")
    return(NULL)
  }

  ind = which(lapply(strsplit(names(g@edgeData),split = "|", fixed = T), function(x){x[[1]][1]}) == from)
  edges = strsplit(names(g@edgeData)[ind], split="|",fixed = T)

  return(edges)
}

edge.exists <- function(g, from, to){
  if(!(from %in% g@nodes)){
    cat("no node ", from, " in the graph\n")
    return(FALSE)
  } else if (!(to %in% g@nodes)){
    cat("no node ", to, " in the graph\n")
    return(FALSE)
  }
  edge.name = paste(from,to,sep = "|")
  ind = which(names(g@edgeData) == edge.name)
  return(length(ind) > 0)
}


get.edge.type <- function(g, from, to){

  edge.name = paste(from,to,sep = "|")
  ind = which(names(g@edgeData) == edge.name)
  if(length(ind) == 0){
    return("none")
  }

  i = ind[1]
  edge = g@edgeData@data[[i]]
  types = c(edge$subtype1, edge$subtype2)
  if(any(grepl("inhibition|repression|dissociation", types)))
    return("inhibition")
  return("activation")
}


