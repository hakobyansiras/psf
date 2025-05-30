
concatenate.summonds <- function(summonds){
  sum = paste(summonds, collapse="+")
}


psf.flow <- function(g, node.ordering, sink.nodes, split = TRUE, sum = FALSE, mult_normalization = FALSE, tmm_mode = FALSE, tmm_updated_mode = FALSE) {

  #   show(i)
  #   k = 0

  node.order <- node.ordering$node.order
  node.rank <- node.ordering$node.rank
  # sink.nodes <- NULL

  nods = names(node.order)
  symb.exprs = vector("list")
  eval.exprs = vector("list")

  #expressions
  E = data.frame(as.numeric(graph::nodeData(g, nods, attr = "expression")), row.names=nods)

  #impacts
  I = matrix(data=NA, nrow=length(nods),ncol= length(nods))
  W = matrix(data = NA, nrow = length(nods), ncol = length(nods))
  rownames(I) = nods
  colnames(I) = nods
  rownames(W) = nods
  colnames(W) = nods

  for (node in nods){
    #     show(node)
    l = length(graph::edgeL(g)[[node]]$edges)
    for (e in 1:l){
      from = node
      to = graph::nodes(g)[graph::edgeL(g)[node][[1]]$edges][e]
      if (!is.na(to) && length(to)>0){
        weight = graph::edgeData(g, from = from, to =to, attr = "weight")[[1]]
        W[from, to] = weight
        impact = graph::edgeData(g, from = from, to =to, attr = "impact")[[1]]
        I[from, to] = impact
      }
    }
  }

  recalc.sinks =F
  if(is.null(sink.nodes)){
    recalc.sinks = T
    cat("\nsink nodes were not supplied. Those will be recalculated ...\n")
  }
  for (i in 1:length(nods)){
    #     k <- k+1
    #     if(names(node.order)[i]=="EntrezGene:10163")
    # show(i)
    parent.nodes <- graph::inEdges(nods[i], g)[[1]]
    
    if(tmm_mode) {
      ## skipping nodes with FC value 1
      parent.nodes <- parent.nodes[which(unlist(graph::nodeData(g, parent.nodes, attr = "signal")) != 1)]
    }
    
    node = nods[i]
    if (length(parent.nodes)>0) {
      #       node.exp <- nodeData(g, i, attr = "expression")[[1]]
      node.exp = E[i,1]
      
      in.signal <- unlist(graph::nodeData(g, parent.nodes, attr = "signal"))
      
      ### Handle parent node signals in loops
      if(any(!(parent.nodes %in% nods[1:i]))) {
        non_processed_parents <- parent.nodes[which(!(parent.nodes %in% nods[1:i]))]
        
        in.signal[which(!(parent.nodes %in% nods[1:i]))] <- unlist(graph::nodeData(g, non_processed_parents, attr = "expression"))
      }
      
      ### input signal processing of function nodes
      if("psf_function" %in% names(graph::nodeDataDefaults(g))) {
        if(unname(unlist(graph::nodeData(g, node, attr = "psf_function"))) %in% c("min", "max", "sum")) {
          in.signal <- get(unname(unlist(graph::nodeData(g, node, attr = "psf_function"))), envir = globalenv())(unlist(graph::nodeData(g, parent.nodes, attr = "signal")))
        }
      }
      
      node.signal <- graph::nodeData(g, nods[i], attr = "signal")[[1]]
      impact = I[parent.nodes,nods[i]]
      weight = W[parent.nodes,nods[i]]
      if(sum){
        node.signal = (node.exp)+sum(in.signal*weight*impact)
      } else {
        if(split){
          proportion = 1
          if(sum(in.signal) != 0){
            proportion = in.signal/sum(in.signal)
            # proportion = 1
            #           node.signal <- (proportion*node.exp)+(in.signal*impact)
            # node.signal <- (node.exp)+sum(in.signal*impact) # no need for proportions in summation - we use this for summinng now
            
            ### previos wroking version
            # node.signal <- sum((proportion*weight*node.exp)*(in.signal^impact)) #esi
            
            ### removes signals from inhibtion parent nodes beacuse of 0^-1 exception
            node.signal <- sum(sapply(1:length(in.signal), function(x) {
              
              if(in.signal[x] == 0 & impact[x] == -1) {
                0
              } else {
                node.exp*(proportion[x]*weight[x]*in.signal[x]^impact[x])
              }
              
            }))
            
            
            # cat("\n node.exp:", node.exp,"node.signal", node.signal, "\n")
            # node.signal <- sum((proportion*node.exp)*(in.signal*impact))
            
            #### special formula for signal exponential decay, not used yet ####
            # a = 2000
            # proportion = in.signal/sum(in.signal)
            # node.signal <- sum((proportion * weight * node.exp) *  a*(2/(1+exp((-2*in.signal^impact)/a)) - 1))
            
          } else {
            #Returns the product of signals without splitting - is for updateing, but applies only in this special case where all the rules are s*t, or 1/s*t
            ##
            ### does this condition similar with self actiavtion? if yes than it should be applied only for those genes which have such a feature.
            ##
            
            proportion = 1
            
            if(all(impact == -1)) {
              node.signal <- node.exp
            } else {
              node.signal <- 0
            }
            
            # node.signal <- sum(node.exp*(weight*in.signal^impact))
            names(node.signal) <- NULL
            
            
            # cat(paste("beging: node",i, "signal:", in.signal, "\n"))
          }
        } else {
          
          if(tmm_mode) {
            if(length(in.signal) > 1) {
              in.signal_adjusted <- weight*in.signal
            } else {
              in.signal_adjusted <- min(weight)*in.signal
              impact <- min(impact)
            }
            
            if(tmm_updated_mode) {
              
              affecting_in_signal <- in.signal_adjusted^impact
              
              node.signal <- prod(c(node.exp, affecting_in_signal))
            } else {
              node.signal <- prod(node.exp*in.signal_adjusted^impact)
            }
            
            
            
          } else {
            #Returns the product of signals without splitting - is for updateing, but applies only in this special case where all the rules are s*t, or 1/s*t
            proportion = 1
            node.signal <- node.exp*weight*prod(in.signal^impact)
          }
          
          
        }
      }
      graph::nodeData(g, nods[i], attr = "signal") <- node.signal

      #       Symb.exprs
      #             in.sign.impacts = vector("list")
      #             signal.base.denoms = vector("list")
      #             for (parent in parent.nodes){
      #               if(!is.null(symb.exprs[[parent]])){
      #                 in.sign.impacts[[parent]] = sprintf("(%s)^(1+I['%s','%s'])", symb.exprs[[parent]], parent, node)
      #                 signal.base.denoms[[parent]] = sprintf('(%s)', symb.exprs[[parent]])
      #               } else{
      #                 in.sign.impacts[[parent]] = sprintf("(E['%s',1])^(1+I['%s','%s'])", parent, parent, node)
      #                 signal.base.denoms[[parent]] = sprintf("(E['%s',1])",parent)
      #               }
      #             }
      #             signal.base.denom = paste(signal.base.denoms, collapse="+")
      #             signal.base = sprintf("E[%d,1]/(%s)", i, signal.base.denom)
      #             in.signal.impact = paste(in.sign.impacts,collapse="+")
      #             symb.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
      #
      #       #Eval.exprs
      in.sign.impacts = vector("list")
      signal.base.denoms = vector("list")
      if(length(parent.nodes) > 1) {
        for (parent in parent.nodes){
          if(!is.null(eval.exprs[[parent]])){
            in.sign.impacts[[parent]] = sprintf("(%s)^(1+I['%s','%s'])", paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
            signal.base.denoms[[parent]] = sprintf("(%s)",paste("eval.exprs[['",parent,"']]",sep=""))
          } else {
            in.sign.impacts[[parent]] = sprintf("(E['%s',1])^(1+I['%s','%s'])", parent, parent, node)
            signal.base.denoms[[parent]] = sprintf("(E['%s',1])",parent)
          }
        }
        signal.base.denom = paste(signal.base.denoms, collapse="+")
        signal.base = sprintf("E[%d,1]/(%s)", i, signal.base.denom)
        in.signal.impact = paste(in.sign.impacts,collapse="+")
        eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
      } else {
        parent = parent.nodes[1]
        if(!is.null(eval.exprs[[parent]])){
          in.signal.impact = sprintf("(%s)^(I['%s','%s'])", paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
        } else {
          in.signal.impact = sprintf("(E['%s',1])^(I['%s','%s'])", parent, parent, node)
        }
        signal.base = sprintf("E[%d,1]", i)
        eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
      }
    } else {
      node.exp = E[i,1]
      graph::nodeData(g, nods[i], attr = "signal") <- node.exp
      #             symb.exprs[[node]] = sprintf("E[%d,1]",i)
      eval.exprs[[node]] = sprintf("E[%d,1]",i)
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
    # sink.nodes <- c(sink.nodes, names(node.order)[i])

  }

  signal.at.sink = NULL
  if(!is.null(sink.nodes)){
    signal.at.sink <- unlist(graph::nodeData(g, sink.nodes, attr = "signal"))
  }
  return(list("graph" = g, "order"=node.ordering, "sink.nodes" = sink.nodes,
              "signal.at.sink" = signal.at.sink,
              "eval.exprs" = eval.exprs, "I" = I, "E" = E))
}


compute.eval.exprs <- function(eval.exprs, E, I, sink.nodes){

  #   start.time = Sys.time()
  for (i in 1:length(eval.exprs)){
    eval.exprs[[i]] = eval(parse(text = eval.exprs[[i]]))
  }
  signal.at.sink = unlist(eval.exprs[sink.nodes])


  #   end.time = Sys.time()
  #   time.eval = end.time - start.time
  return(signal.at.sink)
}


compute.symb.exprs <- function(symb.exprs, E, I, sink.nodes){
  #   start.time = Sys.time()
  signal.at.sink = list()
  for (sink.node in sink.nodes){
    #     show(sink.node)
    signal.at.sink[[sink.node]] = eval(parse(text = symb.exprs[[sink.node]]))
  }


  #   end.time = Sys.time()
  #   time.eval = end.time - start.time
  return(unlist(signal.at.sink))
}


psf.flow.boot <- function(g, FC.matrix, node.ordering, eval.exprs, sink.nodes, I, bst.steps = 200) {
  node.order <- node.ordering$node.order
  node.rank <- node.ordering$node.rank

  genes <- unlist(graph::nodeData(g, attr = "genes"))
  names(genes) <-NULL
  nsamples <- length(genes)
  #   nsamples <- length(nodes(g))
  i = 1

  while (i <=bst.steps){
    #     show(i)
    #     FC.boot <- as.matrix(sample(FC.matrix, nsamples, replace = T))
    #     nodeData(g, attr = "signal") <- FCboot # 1+vector(mode = "numeric", length(FCboot) )
    #     nodeData(g, attr = "expression") <- FCboot
    FC.boot <- matrix(NA, nrow = nsamples, ncol = 1)
    for (j in 1:nrow(FC.boot)) {
      bst.samle <- as.numeric(sample(FC.matrix[j,], 1, replace = T))
      FC.boot[j,] <- bst.samle
    }

    nods = names(node.order)

    rownames(FC.boot) <- genes
    g <- map.gene.data(g, FC.boot)
    E = data.frame(as.numeric(graph::nodeData(g, nods, attr = "expression")), row.names=nods)

    #         boot.signal.at.sink <- compute.symb.exprs(symb.exprs, E, I, sink.nodes)
    boot.signal.at.sink <- compute.eval.exprs(eval.exprs, E, I, sink.nodes)


    if (i == 1) {
      mat.length <- length(boot.signal.at.sink)
      boot.mat <- matrix(NA, nrow = mat.length, ncol = bst.steps)
    }
    boot.mat[,i] <- as.numeric(boot.signal.at.sink)
    i = i+1
  }
  return(boot.mat)
}


sig.calc <- function(signal.at.sink, boot.mat) {
  sig = vector(mode = "numeric", length(signal.at.sink))
  for (i in 1:length(signal.at.sink)) {
    if (signal.at.sink[i] >=1)
      sig[i] = sum(as.numeric(signal.at.sink[i]<=boot.mat[i,]))/length(boot.mat[i,])
    else
      sig[i] = sum(as.numeric(signal.at.sink[i]>=boot.mat[i,]))/length(boot.mat[i,])
  }
  return(sig)
}