psf.flow <- function(g, node.ordering, sink.nodes, split = TRUE, sum = FALSE) {
  
  node.order <- node.ordering$node.order
  node.rank <- node.ordering$node.rank
  # sink.nodes <- NULL
  
  nods = names(node.order)
  symb.exprs = vector("list")
  eval.exprs = vector("list")
  
  #expressions
  E = data.frame(as.numeric(graph::nodeData(g, nods, attr = "expression")), row.names=nods)
  
  #impacts
  I = matrix (data=NA, nrow=length(nods),ncol= length(nods))
  rownames(I) = nods
  colnames(I) = nods
  
  #weights
  W = matrix (data=NA, nrow=length(nods),ncol= length(nods))
  rownames(W) = nods
  colnames(W) = nods
  
  for (node in nods){
    #     show(node)
    l = length(graph::edgeL(g)[[node]]$edges)
    for (e in 1:l){
      from = node
      to = graph::nodes(g)[graph::edgeL(g)[node][[1]]$edges][e]
      if (!is.na(to) && length(to)>0){
        impact = graph::edgeData(g, from = from, to =to, attr = "impact")[[1]]
        weight = graph::edgeData(g, from = from, to =to, attr = "weight.2")[[1]]
        I[from, to] = impact
        W[from, to] = weight
        
      }
    }
  }
  
  recalc.sinks =F
  if(is.null(sink.nodes)){
    recalc.sinks = T
    cat("\nsink nodes were not supplied. Those will be recalculated ...\n")
  }
  
  if(recalc.sinks){
    if (length(parent.nodes) > 0){
      child.nodes <- unlist(graph::edges(g, names(node.order)[i]))
      if (length(child.nodes) >0) {
        child.node.ranks <- node.rank[child.nodes]
        rank.comp = child.node.ranks > node.rank[i]
        if (all(rank.comp))
          sink.nodes <- c(sink.nodes, names(node.order)[i])
      } else {
        sink.nodes <- c(sink.nodes, names(node.order)[i])
      }
    }
  }
  
  
  for (i in 1:length(nods)){
    
    parent.nodes <- graph::inEdges(nods[i], g)[[1]]
    node = nods[i]
    if (length(parent.nodes)>0) {
      #       node.exp <- nodeData(g, i, attr = "expression")[[1]]
      node.exp = E[i,1]
      in.signal <- unlist(graph::nodeData(g, parent.nodes, attr = "signal"))
      
      pi = which(nods %in% parent.nodes)
      #       show(i)
      #       cat(length(parent.nodes), "\n")
      #       in.signal = S[parent.nodes,1]
      node.signal <- graph::nodeData(g, nods[i], attr = "signal")[[1]]
      impact = I[parent.nodes,nods[i]]
      weight = W[parent.nodes,nods[i]]
      if(sum){
        # code when sum is TURE
        # node.signal = (node.exp)+sum(in.signal*impact)
      } else {
        if(split){
          proportion = 1
          if(sum(in.signal) != 0){
            proportion = in.signal/sum(in.signal)
            
            node.signal <- sum((proportion*node.exp)*(weight*in.signal^impact)) #esi
            names(node.signal) <- NULL
          } else {
            #Returns the product of signals without splitting - is for updateing, but applies only in this special case where all the rules are s*t, or 1/s*t
            proportion = 1
            node.signal <- node.exp*(weight*in.signal^impact)
            names(node.signal) <- NULL
            # cat(paste("beging: node",i, "signal:", in.signal, "\n"))
          }
          graph::nodeData(g, nods[i], attr = "signal") <- node.signal
        } else {
          # code where split is FALSE
        }
      }
    } else {
      # assign signal to source nodes
      node.exp = E[i,1]
      graph::nodeData(g, nods[i], attr = "signal") <- node.exp
      #             symb.exprs[[node]] = sprintf("E[%d,1]",i)
      # eval.exprs[[node]] = sprintf("E[%d,1]",i)
    }
    
    
    # sink.nodes <- c(sink.nodes, names(node.order)[i])
    
  }
  
  signal.at.sink = NULL
  if(!is.null(sink.nodes)){
    signal.at.sink <- unlist(graph::nodeData(g, sink.nodes, attr = "signal"))
  }
  return(list("graph" = g, "order"=node.ordering, "sink.nodes" = sink.nodes,
              "signal.at.sink" = signal.at.sink,
              # "eval.exprs" = eval.exprs, 
              "I" = I, "E" = E))
}