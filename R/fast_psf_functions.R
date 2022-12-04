#' Calculates PSF formulas for each node in graphNEL object.
#' @param g graphNEL pathway object.
#' @param node.ordering order of nodes calculated with order.nodes function.
#' @param sink.nodes list of terminal (sink) nodes calculated with determine.sink.nodes function.
#' @param split logical, if true then the incoming signal will be proportionally splitted among the edges.
#' @param sum logical, default value is FALSE. When set to true pathway activity formulas will be caculated via addition, when set to false then activity formulas will be calculated via multiplication.
#' @param tmm_mode when set to true specific PSF configuration will be used for calculation of the pathway activity formulas described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full
#' @param tmm_update_mode when set to true specific PSF configuration will be used for calculation the pathway activity formulas described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full
#' @import graph
eval_formulas <- function(g, node.ordering, sink.nodes, split = TRUE, sum = FALSE, tmm_mode = FALSE, tmm_update_mode = FALSE) {
  
  node.order <- node.ordering$node.order
  node.rank <- node.ordering$node.rank
  
  nods <- names(node.order)
  eval.exprs = vector("list")
  
  #expressions
  # E = data.frame(as.numeric(graph::nodeData(g, nods, attr = "expression")), row.names=nods)
  node_fun = data.frame(as.character(graph::nodeData(g, nods, attr = "psf_function")), row.names=nods)
  
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
    
    parent.nodes <- graph::inEdges(nods[i], g)[[1]]
    
    # if(tmm_mode) {
    #   ## skipping nodes with FC value 1
    #   parent.nodes <- parent.nodes[which(unlist(graph::nodeData(g, parent.nodes, attr = "signal")) != 1)]
    # }
    
    node = nods[i]
    
    if (length(parent.nodes)>0) {
      #       node.exp <- nodeData(g, i, attr = "expression")[[1]]
      # node.exp = E[i,1]
      
      ### input signal processing of function nodes (will be implemented in the future)
      # if("psf_function" %in% names(graph::nodeDataDefaults(g))) {
      #   if(unname(unlist(graph::nodeData(g, node, attr = "psf_function"))) %in% c("min", "max", "sum")) {
      #     in.signal <- get(unname(unlist(graph::nodeData(g, node, attr = "psf_function"))), envir = globalenv())(unlist(graph::nodeData(g, parent.nodes, attr = "signal")))
      #   } else {
      #     in.signal <- unlist(graph::nodeData(g, parent.nodes, attr = "signal"))
      #   }
      # } else {
      #   in.signal <- unlist(graph::nodeData(g, parent.nodes, attr = "signal"))
      # }
      
      if(sum){
        
        in.sign.impacts = vector("list")
        for (parent in parent.nodes){
          if(!is.null(eval.exprs[[parent]])){
            in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(%s)*(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
          } else {
            in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(E['%s',1])*(I['%s','%s'])", parent, node, parent, parent, node)
          }
        }
        in.signal.impact = paste(in.sign.impacts,collapse="+")
        
        eval.exprs[[node]] = sprintf("(E[%d,1])+(%s)", i, in.signal.impact)
        
        
      } else {
        if(split){
          
          in.sign.impacts = vector("list")
          signal.base.denoms = vector("list")
          if(length(parent.nodes) > 1) {
            for (parent in parent.nodes){
              if(!is.null(eval.exprs[[parent]])){
                in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(%s)^(1+I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
                signal.base.denoms[[parent]] = sprintf("(%s)",paste("eval.exprs[['",parent,"']]",sep=""))
              } else {
                in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(E['%s',1])^(1+I['%s','%s'])", parent, node, parent, parent, node)
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
              in.signal.impact = sprintf("(W['%s','%s'])*(%s)^(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
            } else {
              in.signal.impact = sprintf("(W['%s','%s'])*(E['%s',1])^(I['%s','%s'])", parent, node, parent, parent, node)
            }
            signal.base = sprintf("E[%d,1]", i)
            eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
          }
          
          #### special formula for signal exponential decay, not used yet ####
          # a = 2000
          # proportion = in.signal/sum(in.signal)
          # node.signal <- sum((proportion * weight * node.exp) *  a*(2/(1+exp((-2*in.signal^impact)/a)) - 1))
          
          
        } else {

          in.sign.impacts = vector("list")
          signal.base.denoms = vector("list")
          
          if(tmm_mode) {
            
             for (parent in parent.nodes){
               if(!is.null(eval.exprs[[parent]])){
                 in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(%s)^(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
                 signal.base.denoms[[parent]] = sprintf("(%s)",paste("eval.exprs[['",parent,"']]",sep=""))
               } else {
                 in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(E['%s',1])^(I['%s','%s'])", parent, node, parent, parent, node)
                 signal.base.denoms[[parent]] = sprintf("(E['%s',1])",parent)
               }
             }
            
            if(node_fun[i,1] %in% c("min", "max", "sum")) {
              signal.base.denom <- paste(signal.base.denoms, collapse = ", ")
              
              in.signal <- sprintf("%s(%s, %s)", node_fun[i,1], signal.base.denom, "na.rm = TRUE")
              
              eval.exprs[[node]] <- sprintf("(E[%d,1])*(%s)", i, in.signal)
              
            } else {
              
              if(tmm_update_mode) {
                total_signal <- paste(c(sprintf("E[%d,1]", i), in.sign.impacts), collapse = ",")
                
                eval.exprs[[node]] <- sprintf("%s(%s, %s)", "prod", total_signal, "na.rm = TRUE")
                
                # eval.exprs[[node]] <- paste(c(sprintf("E[%d,1]", i), in.sign.impacts), collapse = "*")
              } else {
                total_signal <- paste("(", paste(sprintf("E[%d,1]", i), in.sign.impacts, sep = "*"), ")", collapse = ", ")
                
                eval.exprs[[node]] <- sprintf("%s(%s, %s)", "prod", total_signal, "na.rm = TRUE")
                # eval.exprs[[node]] <- paste("(", paste(sprintf("E[%d,1]", i), in.sign.impacts, sep = "*"), ")", collapse = "*")
                
              }
              
            }
            
            #### on hold, will implemetn this version in the future
            # if(node_fun[i,1] == "sum") {
            # 
            #   for (parent in parent.nodes){
            #     if(!is.null(eval.exprs[[parent]])){
            #       in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(%s)^(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
            #       # signal.base.denoms[[parent]] = sprintf("(%s)",paste("eval.exprs[['",parent,"']]",sep=""))
            #     } else {
            #       in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(E['%s',1])^(I['%s','%s'])", parent, node, parent, parent, node)
            #       # signal.base.denoms[[parent]] = sprintf("(E['%s',1])",parent)
            #     }
            #   }
            # 
            #   in.sign.impact = paste(in.sign.impacts, collapse = "+")
            # 
            # } else {
            #   if(node_fun[i,1] == "max") {
            # 
            #   }
            #   if(node_fun[i,1] == "min")
            # }
            # 
            # 
            # 
            # 
            # sprintf("node_fun['%s', 1]()")
            # 
            # eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
            # 
            # 
            # 
            # if(length(in.signal) > 1) {
            #   in.signal_adjusted <- weight*in.signal
            # } else {
            #   in.signal_adjusted <- min(weight)*in.signal
            #   impact <- min(impact)
            # }
            # 
            # if(tmm_update_mode) {
            # 
            #   affecting_in_signal <- in.signal_adjusted^impact
            # 
            #   node.signal <- prod(c(node.exp, affecting_in_signal))
            # } else {
            #   node.signal <- prod(node.exp*in.signal_adjusted^impact)
            # }
            
            
            
          } else {
            #Returns the product of signals without splitting - is for updateing, but applies only in this special case where all the rules are s*t, or 1/s*t
            
            in.sign.impacts = vector("list")
            if(length(parent.nodes) > 1) {
              for (parent in parent.nodes){
                if(!is.null(eval.exprs[[parent]])){
                  in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(%s)^(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
                } else {
                  in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(E['%s',1])^(I['%s','%s'])", parent, node, parent, parent, node)
                }
              }
              
              in.signal.impact = paste(in.sign.impacts,collapse="*")
              signal.base = sprintf("E[%d,1]", i)
              
              eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
            } else {
              parent = parent.nodes[1]
              if(!is.null(eval.exprs[[parent]])){
                in.signal.impact = sprintf("(W['%s','%s'])*(%s)^(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
              } else {
                in.signal.impact = sprintf("(W['%s','%s'])*(E['%s',1])^(I['%s','%s'])", parent, node, parent, parent, node)
              }
              signal.base = sprintf("E[%d,1]", i)
              eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
            }
            
          }
          
          
        }
      }
      
    } else {
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
  
  return(list("graph" = g, "order"=node.ordering, "sink.nodes" = sink.nodes,
              "eval.exprs" = eval.exprs, "I" = I, "W" = W))
  
}


### use parallel::mclapply function and check if performance is changed, if not use parallel::mclapply

### Implement loop handling in PSF calculation
### LoopDetectR::find_loops(jacobian = igraph::as_adjacency_matrix(igraph::igraph.from.graphNEL(x$graph), sparse = F)) ### loops detection script

#' Calculates pathway activity with PSF algorithm for provided kegg collection based on expression fold change data.
#' @param entrez.fc expression fold change matrix with gene etrez id rownames.
#' @param kegg.collection the list of kegg pathways generated by generate.kegg.collection.from.kgml or generate.kegg.collection function.
#' @param calculate.significance logical, if true then function will also calculate significance for the PSF values by shuffling all the network nodes and checking if the resulted PSF values were calculated by chance.
#' @param bst.steps integer, the number of the interations for shuffling and recalculating PSF values for the significance analysis.
#' @param shuffling_mode charachter, default value is "global". With this argument set of genes is selectited for bootstraping. When set to "global" all frovided FC values will be used to generate randome FC values for the pathway. When set to "local" FC values of pathway genes will be shuffeld each time and used.
#' @param sum logical, default value is FALSE. When set to true pathway activity will be caculated via addition, when set to false then activity willbe calculated via multiplication.
#' @param split logical, if true then the incoming signal will be proportionally splitted among the edges.
#' @param tmm_mode logical, default value is FALSE. When set to true specific PSF configuration will be used for the pathway activity calculation described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full
#' @param tmm_update_mode logical, default value is FALSE. When set to true specific PSF configuration will be used for the pathway activity calculation described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full
#' @param ncores integer, the number of CPUs to use for PSF calculation. Deafult value is 1.
#' @return Kegg collection list with evaluated PSF formulas, calculated PSF values and significance values. Each object of the list corresponds to one pathway where "psf_activities", "p_val_mat", "exp_fc" matrixes are added (each column of the matrix is one sample). This list object can be further used for pathway visualization and PSF report generation.
#' @import graph
#' @import parallel
#' @export
run_psf <- function(entrez.fc, kegg.collection, calculate.significance = F, bst.steps = 200, shuffling_mode = "global", sum = FALSE, split = TRUE, tmm_mode = FALSE, tmm_update_mode = FALSE, ncores = 1) {
  
  ### checking exp fc matrix
  if(!("matrix" %in% class(entrez.fc))) {
    stop(paste0("ERROR: Provided expression data is in ", class(entrez.fc), "class, provide expression matrix in matrix class"))
  }
  
  if(any(is.na(entrez.fc))) {
    stop("ERROR: expression matrix contain NA values")
  }
  
  if(any(entrez.fc == 0)) {
    stop("ERROR: expression matrix contain 0 values, make sure to provide expression fold change values")
  }
  
  
  ### runing psf on pathway
  psf_processing <- function(entrez.fc, pathway) {
    
    ### building psf formulas
    cat(paste0("Evaluating PSF formulas for ", pathway$attrs$title, "\n"))
    
    eval_g <- eval_formulas(g = pathway$graph, node.ordering = pathway$order, sink.nodes = pathway$sink.nodes, split = split, sum = sum, tmm_mode = tmm_mode, tmm_update_mode = tmm_update_mode)
    
    cat(paste0("Processing Expression FC data for ", pathway$attrs$title, "\n"))
    
    ### building exp FC matrix for the pathway
    gene.data <- graph::nodeData(pathway$graph)[names(pathway$order$node.order)]
    pathway_exp_data <- lapply(gene.data, function(x) {
      
      genes.in.node = which(rownames(entrez.fc) %in% x$genes)
      expression.values <- entrez.fc[genes.in.node,, drop = F]
      
      if (nrow(expression.values)>0){
        if (x$type == "gene"){
          expression <- matrix(colMeans(expression.values), ncol = ncol(entrez.fc))
        } else { 
          if(x$type == "group") {
            expression <- apply(expression.values, 2, FUN = min)
          } else {
            expression <- matrix(rep(1, ncol(entrez.fc)), ncol = ncol(entrez.fc))
          }
        }
      } else {
        expressions <- matrix(rep(1, ncol(entrez.fc)), ncol = ncol(entrez.fc))
      }
      
    })
    
    pathway_exp_data <- Reduce(rbind, pathway_exp_data)
    
    rownames(pathway_exp_data) <- names(pathway$order$node.order)
    colnames(pathway_exp_data) <- colnames(entrez.fc)
    
    
    ### calculating psf with pathway formulas
    I = eval_g$I
    W = eval_g$W
    eval.exprs <- eval_g$eval.exprs
    
    eval_calculator <- function(exp_ind, exp_mat) {
      E <- exp_mat[,exp_ind, drop = F]
      for (i in 1:length(eval.exprs)){
        eval.exprs[[i]] = eval(parse(text = eval.exprs[[i]]))
      }
      return(eval.exprs)
    }
    
    psf_activities <- simplify2array(parallel::mclapply(1:ncol(pathway_exp_data), function(x) {eval_calculator(exp_ind = x, exp_mat = pathway_exp_data)}, mc.cores = ncores))
    
    # single threaded version
    # psf_activities <- sapply(1:ncol(pathway_exp_data), function(x) {
    #   eval_calculator(x)
    # })
    
    ### not needed for bare psf matrix output
    # node_labels <- unlist(graph::nodeData(pathway$graph, rownames(pathway_exp_data), attr = "label"))
    # node_ids <- rownames(pathway_exp_data)
    
    ### calculating significance (need to be iproved, very slow with this setup)
    if(calculate.significance) {
      p_val_mat <- Reduce(cbind, lapply(1:ncol(pathway_exp_data), function(x) {
        
        if(shuffling_mode == "local") {
          random_fc_mat <- sapply(1:bst.steps, function(y) {
            pathway_exp_data[sample(1:nrow(pathway_exp_data), size = nrow(pathway_exp_data)),x, drop = F]
          })
          rownames(random_fc_mat) <- rownames(pathway_exp_data)
        } else {
          random_fc_mat <- sapply(1:bst.steps, function(y) {
            entrez.fc[sample(1:nrow(entrez.fc), size = nrow(pathway_exp_data)),x, drop = F]
          })
          rownames(random_fc_mat) <- rownames(pathway_exp_data)
        }
        
        
        psf_activities_dist <- simplify2array(parallel::mclapply(1:ncol(random_fc_mat), function(x) {eval_calculator(exp_ind = x, exp_mat = random_fc_mat)}, mc.cores = ncores))
        
        
        sapply(rownames(psf_activities_dist), function(z) {
          if(psf_activities[z,x] >= 1) {
            sum(unlist(psf_activities_dist[z,]) > unlist(psf_activities[z,x]))/bst.steps
          } else {
            sum(unlist(psf_activities_dist[z,]) < unlist(psf_activities[z,x]))/bst.steps
          }
        })
        
      }))
      
      psf_activities <- apply(psf_activities, 2, as.numeric)
      colnames(psf_activities) <- colnames(pathway_exp_data)
      rownames(psf_activities) <- rownames(pathway_exp_data)
      
      p_val_mat <- apply(as.matrix(p_val_mat), 2, as.numeric)
      colnames(p_val_mat) <- colnames(pathway_exp_data)
      rownames(p_val_mat) <- rownames(pathway_exp_data)
      
      ### not used, bare matriexes are preffered for downstream anaysis
      # psf_and_pvals <- Reduce(cbind, lapply(1:ncol(psf_activities), function(x) {
      #   cbind(psf_activities[,x], as.matrix(p_val_mat)[,x])
      # }))
      # 
      # psf_and_pvals <- apply(psf_and_pvals, 2, as.numeric)
      # 
      # colnames(psf_and_pvals) <- unlist(lapply(colnames(pathway_exp_data), function(x) {c(x, paste0(x, "_pval"))}))
      # 
      # psf_and_pvals <- data.frame(node_ids, node_labels, psf_and_pvals, stringsAsFactors = F, row.names = node_ids)
      # 
      # pathway$psf_and_pvals <- psf_and_pvals
      
      pathway$psf_activities <- psf_activities
      
      pathway$p_val_mat <- p_val_mat
      
      pathway$exp_fc <- pathway_exp_data
      
    } else {
      psf_activities <- apply(psf_activities, 2, as.numeric)
      
      colnames(psf_activities) <- colnames(pathway_exp_data)
      
      rownames(psf_activities) <- rownames(pathway_exp_data)
      # psf_activities <- data.frame(node_ids, node_labels, psf_activities, stringsAsFactors = F, row.names = node_ids)
      
      pathway$psf_activities <- psf_activities
      
      pathway$exp_fc <- pathway_exp_data
    }
    
    
    pathway$eval.exprs <- eval_g$eval.exprs
    
    pathway$I <- I
    pathway$W <- W
    
    return(pathway)
    
  }
  
  return(lapply(kegg.collection, function(x) {
    psf_processing(entrez.fc, x)
  }))
  
}



### complete implementation of specific influenced node(s) calculation (issues, single influence node, many samples with differnet influence node order, etc.)
#' Performs partial influence analysis which evaluates effect of each pathway node on specific node(s) of the pathway.
#' @param pathway Pathway onbject with calculated PSF activities with run_psf function. 
#' @param influenced_node sigle id of node or vector of node ids based on which partial influencse will be calculated.
#' @param influence_direction in which direction node affects the signal of target node(s). Possible values c("+", "-", "any"). Default value is "any".
#' @param node_combinations number of node combinations for partial influence analysis. Note: number of PSF calculations increases by the exponent of combinations(N nodes^combinations). Default value is 1.
#' @param get_influence_matrix When set to true influence matrix will be returned instead of influencial nodes. Each column in influence matrix represent the log ratio of default and neutralized node(s) psf profile.
#' @param ncores integer, the number of CPUs to use for calculation. Deafult value is 1.
#' @return Returns list with 2 objects. "ordered_influence_nodes" matrix with ordered nodes based on their influence where each column correspond to one sample. "mean_influence_params" dataframe with mean infulence parameters for each node (NULL when analysis performed for 1 sample).
#' @import parallel
#' @import graph
#' @export
run_pi <- function(pathway, influenced_node, influence_direction = "any", node_combinations = 1, get_influence_matrix = FALSE, ncores = 1) {
  
  if(!("eval.exprs" %in% names(pathway))) {
    stop("Please provide pathway with evaluated activity")
  }
  
  node_combs <- combn(graph::nodes(pathway$graph), node_combinations, simplify = F)
  
  eval.exprs <- pathway$eval.exprs
  I <- pathway$I
  W <- pathway$W
  
  eval_calculator <- function(exp_ind, exp_data) {
    E <- exp_data[,exp_ind, drop = F]
    for (i in 1:length(eval.exprs)){
      eval.exprs[[i]] = eval(parse(text = eval.exprs[[i]]))
    }
    return(eval.exprs)
  }
  
  psf_iterator <- function(node_combs, exp_mat, exp_ind) {
    
    updated_exp <- exp_mat
    
    updated_exp[node_combs,] <- 1
    
    return(eval_calculator(exp_ind = exp_ind, exp_data = updated_exp))
  }
  
  
  psf_influence_list <- lapply(1:ncol(pathway$exp_fc), function(x) {
    
    influence_psf_mat <- Reduce(cbind, parallel::mclapply(node_combs, function(y) {psf_iterator(node_combs = y, exp_mat = pathway$exp_fc, exp_ind = x)}, mc.cores = ncores))
    
    influence_psf_mat <- apply(influence_psf_mat, 2, as.numeric)
    
    colnames(influence_psf_mat) <- sapply(node_combs, function(y) {paste(y, collapse = "_")})
    rownames(influence_psf_mat) <- rownames(pathway$exp_fc)
    
    psf_difference <- log(influence_psf_mat/pathway$psf_activities[,x])
    
    if(get_influence_matrix) {
      psf_difference
    } else {
      if(influence_direction == "any") {
        influence_nodes_order <- order(abs(colSums(psf_difference[influenced_node,, drop = FALSE]))[which(abs(colSums(psf_difference[influenced_node,, drop = FALSE])) > 0)], decreasing = T)
        
        ordered_influence_nodes <- colSums(psf_difference[influenced_node,, drop = FALSE])[which(abs(colSums(psf_difference[influenced_node,, drop = FALSE])) > 0)][influence_nodes_order]
        
        # ordered_influence_nodes <- sort(abs(colSums(psf_difference[influenced_node,, drop = FALSE]))[which(abs(colSums(psf_difference[influenced_node,, drop = FALSE])) > 0)], decreasing = T)
      }
      
      if(influence_direction == "+") {
        ordered_influence_nodes <- sort(colSums(psf_difference[influenced_node,, drop = FALSE])[which(colSums(psf_difference[influenced_node,, drop = FALSE]) > 0)], decreasing = T)
      }
      
      if(influence_direction == "-") {
        ordered_influence_nodes <- sort(colSums(psf_difference[influenced_node,, drop = FALSE])[which(colSums(psf_difference[influenced_node,, drop = FALSE]) < 0)])
      }
      
      influence_nodes_dataframe <- data.frame(node_id = names(ordered_influence_nodes), inf_ratio = unname(ordered_influence_nodes), stringsAsFactors = F)
      colnames(influence_nodes_dataframe) <- paste0(colnames(pathway$exp_fc)[x], "_", colnames(influence_nodes_dataframe))
      influence_nodes_dataframe
    }
  })
  
  if(get_influence_matrix) {
    names(psf_influence_list) <- colnames(pathway$exp_fc)
    return(psf_influence_list)
  } else {
    
    
    if(ncol(pathway$exp_fc) > 1) {
      mutual_inf_nodes <- Reduce(intersect, lapply(psf_influence_list, function(y) {y[,1]}))
      
      mean_node_rank <- sort(sapply(mutual_inf_nodes, function(x) {
        mean(sapply(psf_influence_list, function(y) {which(y[,1] == x)}))
      }))
      
      mean_node_influence_ratio <- sapply(mutual_inf_nodes, function(x) {
        mean(sapply(psf_influence_list, function(y) {y[which(y[,1] == x), 2]}))
      })
      
      mean_influence_params <- data.frame(node_id = names(mean_node_rank), 
                                          mean_node_rank = mean_node_rank, 
                                          mean_inf_ratio = mean_node_influence_ratio[names(mean_node_rank)], 
                                          stringsAsFactors = F)
    }
    
    ordered_influence_nodes <- Reduce(cbind, psf_influence_list)
    
    return(list(ordered_influence_nodes = ordered_influence_nodes, mean_influence_params = mean_influence_params))
  }
  
}


#' Plots the pathway with colored nodes and labels
#' @param pathway pathway object
#' @param plot_type network visualization type. Possible values c("kegg", "visnet"). When "kegg" option is used the pathway will be plotted over kegg png image. With "visnet" option, function will plot interactive network.
#' @param color_nodes type of node values to be visualized. When value type is specified pathway nodes will be color coded with expression FC values or PSF values and color legend will be added to the pathway plot. Possible values are c(NULL, "psf_activities", "exp_fc"). Default value is NULL.
#' @param sample_id name of the sample which will be used to visualize psf or exp fc values on pathway. To get averaged values accross samples set value to "mean". Default value is "mean". 
#' @param log_norm log transform PSF and expression values before color mapping. Default value is TRUE.
#' @param use_old_images use old kegg images(for use with curated pathway collection). Default value is FALSE.
#' @param highlight_nodes single value of node id or a vector of ids to be highlighted in plotted pathway. Default value is NULL.
#' @param highlight_color Highlighted nodes color(s). Default values is "red".
#' @param y_adj_text Numeric value to adjust Y position of node labels. Depending on R version this value should be adjusted to have accurate text alignment. Default value is 0.
#' @param y_adj_sink Numeric value to adjust Y position of sink node sign labels. Depending on R version this value should be adjusted to have accurate alignment. Default value is 0.
#' @import graph
#' @import visNetwork
#' @import RCurl
#' @export
plot_pathway <- function(pathway, plot_type = "visnet",
                         color_nodes = NULL, sample_id = "mean", log_norm = T, 
                         highlight_nodes = NULL, highlight_color = "red", # hover_data = "all"
                         y_adj_text = 0, y_adj_sink = 0, use_old_images = F) {
  
  ### chechikng highlight nodes and their border colors
  if(length(highlight_color) > 1) {
    if(length(highlight_color) == length(highlight_nodes)) {
      highlight_color_vector <- setNames(object = highlight_color, nm = highlight_nodes)
    } else {
      stop("Error: highlighted nodes and their colors must be in the same length")
    }
  } else {
    highlight_color_vector <- setNames(object = rep(highlight_color, length(highlight_nodes)), nm = highlight_nodes)
  }
  
  
  if("exp_fc" %in% names(pathway)) {
    mapping_data = TRUE
    if(sample_id == "mean") {
      pathway_exp_values <- rowMeans(pathway$exp_fc)
      pathway_psf_values <- rowMeans(pathway$psf_activities)
    } else {
      pathway_exp_values <- pathway$exp_fc[,sample_id]
      pathway_psf_values <- pathway$psf_activities[,sample_id]
    }
    
    if(log_norm) {
      pathway_exp_values <- round(log(drop(pathway_exp_values)[order(drop(pathway_exp_values))] + 0.00001), digits = 5)
      pathway_psf_values <- round(log(drop(pathway_psf_values)[order(drop(pathway_psf_values))] + 0.00001), digits = 5)
    } else {
      pathway_exp_values <- round(drop(pathway_exp_values)[order(drop(pathway_exp_values))], digits = 5)
      pathway_psf_values <- round(drop(pathway_psf_values)[order(drop(pathway_psf_values))], digits = 5)
    }
  } else {
    mapping_data = FALSE
  }
  
  if(!is.null(color_nodes)) {
    
    if(mapping_data) {
      if(color_nodes == "psf_activities") {
        pathway_node_values <- pathway_psf_values
      } else {
        pathway_node_values <- pathway_exp_values
      }
      
      node_colors <- color_code(values = pathway_node_values, pal1 = pal1, pal2 = pal2, log_scale = log_norm)
      
      node_colors <- data.frame(node_id = names(pathway_node_values)[c(which(pathway_node_values <= 0), which(pathway_node_values > 0))], 
                                col = node_colors,
                                text_col = unname(sapply(node_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                                stringsAsFactors = F
      )
      
      node_colors <- node_colors[which(node_colors$node_id %in% graph::nodes(pathway$graph)[which(unlist(graph::nodeData(pathway$graph, attr = "type")) != "map")]),]
      
    } else {
      stop("Please provide pathway with evalueated activity")
    }  
    
  } else {
    node_colors <- NULL
  }
  
  
  graphical_data <- graphical_data_generator(pathway)
  
  node_graphics <- graphical_data$node_coords
  
  
  if(color_nodes == "psf_activities") {
    col_legend_title = ifelse(log_norm, "Log PSF value", "PSF value")
  } else {
    col_legend_title = ifelse(log_norm, "Log FC value", "FC value")
  }
  
  if(plot_type == "kegg") {
    if(use_old_images) {
      img_path <- system.file("extdata", "old_imgs", paste0(gsub("path:", "", pathway$attrs$name), ".png"), package="psf")
    } else {
      if(file.exists(system.file("extdata", "pathway_imgs", paste0(gsub("path:", "", pathway$attrs$name), ".png"), package="psf"))) {
        img_path <- system.file("extdata", "pathway_imgs", paste0(gsub("path:", "", pathway$attrs$name), ".png"), package="psf")
      } else {
        # image <- paste0("http://rest.kegg.jp/get/", pathway_id, "/image")
        
        img_path <- tempfile()
        download.file(pathway$attrs$image, destfile = img_path, method = "curl")
      }
      
    }
    
    pathway_image = magick::image_read(img_path)
    
    img <- magick::image_draw(pathway_image)
    
    
    ### node exp coloring
    node_graphics$x_center <- node_graphics$x_start + (node_graphics$x_end - node_graphics$x_start)/2
    
    node_graphics$y_center <- node_graphics$y_start + (node_graphics$y_end - node_graphics$y_start)/2
    
    rownames(node_colors) <- node_colors$node_id
    
    rownames(node_graphics) <- node_graphics$node_id
    
    sink_node_graphics <- node_graphics[which(node_graphics$sink),]
    
    
    if(any(node_graphics$node_id %in% node_colors$node_id)) {
      coloring_set <- node_graphics[node_colors$node_id,]
      graphics::rect(coloring_set$x_start,
                     coloring_set$y_start - 1,
                     coloring_set$x_end,
                     coloring_set$y_end,
                     # border = coloring_set$border_color,
                     border = NA,
                     # lty = coloring_set$lty_type,
                     lwd=2,
                     col = grDevices::adjustcolor(node_colors$col, alpha.f = 1)
      )
      
      graphics::text(x = coloring_set$x_center,
                     y = coloring_set$y_start + y_adj_text,
                     labels = coloring_set$gr_name,
                     col = node_colors$text_col, adj = c(0,0.2) + c(0.48, 1))
      
    }
    
    ### adding sink node labels
    graphics::text(x = sink_node_graphics$x_end + 10,
                   y = sink_node_graphics$y_center - 30 + y_adj_sink, cex = 3,
                   labels = rep("*", nrow(sink_node_graphics)),
                   col = rep("#9ACD32", nrow(sink_node_graphics)), adj = c(0,0.2) + c(0.48, 1))
    
    if(!is.null(highlight_nodes)) {
      highlight_set <- node_graphics[which(node_graphics[,"node_id"] %in% highlight_nodes),]
      
      rect( highlight_set$x_start, 
            highlight_set$y_start, 
            highlight_set$x_end, 
            highlight_set$y_end, 
            border = highlight_color_vector[highlight_set$node_id], lty = "solid", lwd=2, col = adjustcolor( "#a3297a", alpha.f = 0))
    }
    
    
    
    
    ### scale color bar
    if(!is.null(color_nodes)) {
      color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = range(pathway_node_values), digits=3, prompt=FALSE,
                         lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    }
    
    text(x = c(magick::image_info(img)$width - 88, magick::image_info(img)$width - 30),
         y = c(70, 65 + y_adj_sink - 6), cex = c(1.5, 3),
         labels = c("Sink node", "*"),
         col = c("#000000", "#9ACD32"), adj = c(0,0.2) + c(0.48, 1))
    
    
    ## color grop nodes
    if(length(pathway$group_nodes) > 0 ) {
      lapply(pathway$group_nodes, function(z) {
        graphics::rect( z$kegg.gr.x-z$kegg.gr.width*0.5, 
                        z$kegg.gr.y+z$kegg.gr.height*0.5, 
                        z$kegg.gr.x+z$kegg.gr.width*0.5, 
                        z$kegg.gr.y-z$kegg.gr.height*0.5, 
                        border = "yellow", lty = "dashed", lwd=2)
      })
    }
    
    
    dev.off()
    
    return(img)
  } else {
    
    graphical_data$edge_coords <- graphical_data$edge_coords[which(graphical_data$edge_coords$lty == "solid"),]
    graphical_data$node_coords <- graphical_data$node_coords[which(graphical_data$node_coords$exist),]
    
    node_shapes <- unname(sapply(graphical_data$node_coords$node_name, function(x) {
      if(grepl("cpd",x)) {
        "dot"
      } else {
        "box"
      }
    }))
    
    if(!is.null(node_colors)) {
      color <- unname(sapply(graphical_data$node_coords$node_id, function(x) {
        if(x %in% node_colors$node_id) {
          node_colors[which(node_colors$node_id == x),"col"]
        } else {
          "#BFFFBF"
        }
      }))
      
      font_color <- unname(sapply(graphical_data$node_coords$node_id, function(x) {
        if(x %in% node_colors$node_id) {
          node_colors[which(node_colors$node_id == x),"text_col"]
        } else {
          "#000000"
        }
      }))
      
      if(log_norm) {
        hover_string <- paste(graphical_data$node_coords$hover_name, 
                              paste("Node id", graphical_data$node_coords$node_id),
                              paste("Log exp FC", pathway_exp_values[graphical_data$node_coords$node_id]),
                              paste("Log PSF", pathway_psf_values[graphical_data$node_coords$node_id]),
                              sep = "<br>")
      } else {
        hover_string <- paste(graphical_data$node_coords$hover_name, 
                              paste("Node id", graphical_data$node_coords$node_id),
                              paste("Exp FC", pathway_exp_values[graphical_data$node_coords$node_id]),
                              paste("PSF", pathway_psf_values[graphical_data$node_coords$node_id]),
                              sep = "<br>")
      }
      
    } else {
      
      if(mapping_data) {
        hover_string <- paste(graphical_data$node_coords$hover_name, 
                              paste("Node id", graphical_data$node_coords$node_id),
                              paste("Exp FC", pathway_exp_values[graphical_data$node_coords$node_id]),
                              paste("PSF", pathway_psf_values[graphical_data$node_coords$node_id]),
                              sep = "<br>")
      } else {
        hover_string <- paste(graphical_data$node_coords$hover_name, 
                              paste("Node id", graphical_data$node_coords$node_id),
                              sep = "<br>")
      }
      
      
      
      color <- rep("#BFFFBF", nrow(graphical_data$node_coords))
      
      font_color <- rep("#000000", nrow(graphical_data$node_coords))
      
    }
    
    if(is.null(highlight_nodes)) {
      border_color <- unname(sapply(graphical_data$node_coords$sink, function(x) {
        if(x) {
          "#0099cc"
        } else {
          "#BFFFBF"
        }
      }))
    } else {
      
      border_color <- unname(sapply(graphical_data$node_coords$node_id, function(x) {
        if(x %in% highlight_nodes) {
          if(length(highlight_color) > 1) {
            highlight_color_vector[x]
          } else {
            highlight_color
          }
        } else {
          "#BFFFBF"
        }
      }))
      
    }
    
    size <- unname(sapply(graphical_data$node_coords$node_name, function(x) {
      if(grepl("cpd",x)) {
        10
      } else {
        25
      }
    }))
    
    nodes <- data.frame(id = graphical_data$node_coords$node_id,
                        image = rep("unselected", nrow(graphical_data$node_coords)),
                        label = graphical_data$node_coords$gr_name,
                        shape = node_shapes, color.background = color, 
                        color.border = border_color, 
                        # color.highlight = color,
                        borderWidth = 2,
                        title = hover_string,
                        font.size = rep(22, nrow(graphical_data$node_coords)), size = size,
                        font.color = font_color,
                        x = (graphical_data$node_coords$x_start + graphical_data$node_coords$x_end)/2,
                        y = (graphical_data$node_coords$y_start + graphical_data$node_coords$y_end)/2
    )
    
    if(!is.null(node_colors)) {
      
      magick::autoviewer_disable() ### to avoid legend plotting before netowrk rendering
      legend_img <- magick::image_device(width = 480, height = 480)
      plot.new()
      
      color_legend_maker(x = 0.05, y = 0, leg = 0.9, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = range(pathway_node_values), digits=3, prompt=FALSE,
                         lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
      
      temp_legend <- tempfile()
      
      magick::image_write(magick::image_trim(legend_img, fuzz = 0), path = temp_legend)
      
      legend_path <- paste('data:image/png;base64', RCurl::base64Encode(readBin(temp_legend, 'raw', file.info(temp_legend)[1, 'size']), 'txt'), sep = ',')
      
      legend_data_frame <- data.frame(
        id = as.character(max(as.integer(graphical_data$node_coords$node_id)) + 1),
        image = legend_path,
        label = "Color legend", shape = "image", color.background = "",
        color.border = "", borderWidth = 0, title = "",
        font.size = 32, size = 30, font.color = "#000000",
        x = max(graphical_data$node_coords$x_end),
        y = min(graphical_data$node_coords$y_end) - 10
      )
      
      nodes <- rbind(nodes, legend_data_frame)
    }
    
    arrows_type <- c("arrow","bar")
    names(arrows_type) <- c("simple", "T")
    
    edges <- data.frame(from = graphical_data$edge_coords$from, to = graphical_data$edge_coords$to,
                        color = graphical_data$edge_coords$col,
                        arrows.to.enabled = rep(TRUE, length(graphical_data$edge_coords$col)),
                        arrows.to.type = arrows_type[graphical_data$edge_coords$arr.type]
    )
    
    saving_pathway_name <- gsub("(", "", pathway$attrs$title, fixed = T)
    saving_pathway_name <- gsub(")", "", saving_pathway_name, fixed = T)
    saving_pathway_name <- gsub("-", "_", saving_pathway_name)
    saving_pathway_name <- gsub(" ", "_", saving_pathway_name)
    saving_pathway_name <- gsub("___", "_", saving_pathway_name)
    
    visNetwork(nodes = nodes, edges = edges, width = "100%", height = "800px") %>% 
      visIgraphLayout(layout = "layout_nicely") %>% 
      visExport(name = saving_pathway_name) %>%
      visInteraction(navigationButtons = TRUE, multiselect = T)
    
  }
  
}


### report genrating function
generate_psf_report <- function(psf_list, folder_name) {
  
  test_psf_result
  
  lapply(psf_list, function(x) {
    
    plot_pathway
    
  })
  
}

#' Converts graphNEL object to 2 data frames(node_table, edge_table)
#' @param pathway pathway list object generated with package.
#' @param extended set to TRUE to get more detailed data frames with graphical information. Default value is FALSE.
#' @import graph
#' @export
graphnel_to_df <- function(pathway, extended = FALSE) {
  
  splitted_interactions <- strsplit(names(pathway$graph@edgeData@data), split = "|", fixed = T)
  from <- as.character(sapply(splitted_interactions, "[[", 1))
  to <- as.character(sapply(splitted_interactions, "[[", 2))
  
  if(extended) {
    entrez_id <- unname(sapply(pathway$graph@nodes, function(y) {
      ifelse(is.null(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),
             as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))), 
             paste(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), collapse = ","))
    }))
    
    hover_name <- unname(sapply(pathway$graph@nodes, function(y) {
      ifelse(is.null(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),
             paste0(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))), " (",
                    unlist(strsplit(kegg_compounds_to_full_name[as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))),], split = ";"))[1],
                    ")"
             ),
             paste0(
               paste0(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), 
                      paste0("(", entrez_to_symbol[as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),], ")")),
               collapse = ", "
             )
      )
    }))
    
    node_shapes <- unname(sapply(as.character(unlist(graph::nodeData(pathway$graph, attr = "kegg.name"))), function(x) {
      if(grepl("cpd",x)) {
        "dot"
      } else {
        "box"
      }
    }))
    
    
    color <- unname(sapply((unlist(graph::nodeData(pathway$graph, attr = "kegg.id")) %in% pathway$sink.nodes), function(x) {
      if(x) {
        "#0099cc"
      } else {
        "#BFFFBF"
      }
    }))
    
    size <- unname(sapply(as.character(unlist(graph::nodeData(pathway$graph, attr = "kegg.name"))), function(x) {
      if(grepl("cpd",x)) {
        10
      } else {
        25
      }
    }))
    
    sink = (unlist(graph::nodeData(pathway$graph, attr = "kegg.id")) %in% pathway$sink.nodes)
    
    node_table <- data.frame(id = as.character(unlist(graph::nodeData(pathway$graph, attr = "kegg.id"))),
                             label = as.character(unlist(graph::nodeData(pathway$graph, attr = "label"))),
                             component_id_s = unname(sapply(pathway$graph@nodes, function(y) {
                               ifelse(is.null(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),
                                      as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))), 
                                      paste(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), collapse = ","))
                             })),
                             type = as.character(unlist(graph::nodeData(pathway$graph, attr = "type"))),
                             psf_function = as.character(unlist(graph::nodeData(pathway$graph, attr = "psf_function"))),
                             exp_fc = as.numeric(unlist(graph::nodeData(pathway$graph, attr = "expression"))),
                             signal = as.numeric(unlist(graph::nodeData(pathway$graph, attr = "signal"))),
                             title = hover_name,
                             shape = node_shapes,
                             color = color,
                             color.border = "#BFFFBF",
                             borderWidth = 2,
                             font.size = 22,
                             size = size,
                             font.color = "#000000",
                             x = as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.x"))),
                             y = as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.y"))),
                             row.names = as.character(unlist(graph::nodeData(pathway$graph, attr = "kegg.id"))),
                             stringsAsFactors = F
    )
    
    kegg_arrows_type <- c("arrow", "arrow", "arrow", "arrow", "bar", "arrow", "arrow", "arrow", "bar", "arrow", "", "arrow", "arrow", "bar", "arrow", "arrow")
    names(kegg_arrows_type) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
    
    line_col <- c("red", "red", "red", "red", "blue","red", "red", "red", "blue", "red", "red", "red", "red", "blue", "red", "red")
    names(line_col) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
    
    edge_table <- data.frame(from = from,
                             to = to,
                             color = line_col[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "subtype1")))],
                             arrows.to.enabled = TRUE,
                             arrows.to.type = kegg_arrows_type[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "subtype1")))],
                             type = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype1")),
                             subtype = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype2")),
                             state = "",
                             weight = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "weight")),
                             stringsAsFactors = F
    )
    
  } else {
    
    node_table <- data.frame(node_id = as.character(unlist(graph::nodeData(pathway$graph, attr = "kegg.id"))),
                             label = as.character(unlist(graph::nodeData(pathway$graph, attr = "label"))),
                             component_id_s = unname(sapply(pathway$graph@nodes, function(y) {
                               ifelse(is.null(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),
                                      as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))), 
                                      paste(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), collapse = ","))
                             })),
                             type = as.character(unlist(graph::nodeData(pathway$graph, attr = "type"))),
                             psf_function = as.character(unlist(graph::nodeData(pathway$graph, attr = "psf_function"))),
                             exp_fc = as.numeric(unlist(graph::nodeData(pathway$graph, attr = "expression"))),
                             signal = as.numeric(unlist(graph::nodeData(pathway$graph, attr = "signal"))),
                             x = as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.x"))),
                             y = as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.y"))),
                             stringsAsFactors = F
    )
    
    edge_table <- data.frame(from = from,
                             to = to,
                             type = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype1")),
                             subtype = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype2")),
                             state = "",
                             weight = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "weight")),
                             stringsAsFactors = F
    )
    
  }
  
  return(list(node_table = node_table, edge_table = edge_table))
  
}

#' Converts 2 data frames(node_table, edge_table) to graphNEL object
#' @param node_table node data frame with specified columns. Mandatory columns: label, component_id_s, type. Non mandatory columns: psf_function, x, y.
#' @param edge_table edge data frame with specified columns. Mandatory columns: from, to, type. Non mandatory columns subtype, state, weight.
#' @importFrom "igraph" "as_graphnel"
#' @importFrom "igraph" "graph_from_data_frame"
#' @import graph
#' @export
df_to_graphnel <- function(node_table, edge_table) {
  
  # creating graphNEL object
  g <- as_graphnel(
    graph_from_data_frame(edge_table[,c("from", "to")], directed = T, vertices = node_table$node_id)
  )
  
  ### setting node default values
  graph::nodeDataDefaults(g,attr="genes") <- 0
  graph::nodeDataDefaults(g,attr="expression") <- 1
  graph::nodeDataDefaults(g,attr="signal") <- 1
  graph::nodeDataDefaults(g,attr="type") <- "gene"
  graph::nodeDataDefaults(g,attr="label") <- ""
  graph::nodeDataDefaults(g,attr="components") <- "NA"
  graph::nodeDataDefaults(g,attr="psf_function") <- "mean"
  graph::nodeDataDefaults(g,attr="kegg.gr.x") <- 0
  graph::nodeDataDefaults(g,attr="kegg.gr.y") <- 0
  graph::nodeDataDefaults(g,attr="kegg.id") <- "0"
  
  ### filling node values
  graph::nodeData(g, attr = "genes") <- node_table$component_id_s
  graph::nodeData(g, attr = "label") <- node_table$label
  graph::nodeData(g, attr = "type") <- node_table$type
  graph::nodeData(g, attr = "psf_function") <- node_table$psf_function
  graph::nodeData(g, attr = "kegg.gr.x") <- node_table$x
  graph::nodeData(g, attr = "kegg.gr.y") <- node_table$y
  
  
  edge.attrs = list(impact="impact", type="type", subtype1="subtype1", subtype2="subtype2",  weight = "weight")
  
  graph::edgeDataDefaults(g,attr=edge.attrs$impact) <- 1
  graph::edgeDataDefaults(g,attr=edge.attrs$weight) <- edge_table$weight
  
  graph::edgeDataDefaults(g,attr=edge.attrs$type) <- edge_table$type
  graph::edgeDataDefaults(g,attr=edge.attrs$subtype1) <- edge_table$subtype
  graph::edgeDataDefaults(g,attr=edge.attrs$subtype2) <- ""
  
  pathway <- list(graph = g)
  pathway$sink.nodes <- psf::determine.sink.nodes(pathway)
  pathway$order <-  psf::order.nodes(pathway$graph)
  pathway$graph <- psf::set.edge.impacts(pathway$graph)
  return(pathway)
  
}
