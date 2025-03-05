#' Calculates PSF formulas for each node in graphNEL object.
#' @param g graphNEL pathway object.
#' @param node.ordering order of nodes calculated with order.nodes function.
#' @param sink.nodes list of terminal (sink) nodes calculated with determine.sink.nodes function.
#' @param split logical, if true then the incoming signal will be proportionally splitted among the edges.
#' @param sum logical, default value is FALSE. When set to true pathway activity formulas will be calculated via addition, when set to false then activity formulas will be calculated via multiplication.
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
#' @param bst.steps integer, the number of the iterations for shuffling and recalculating PSF values for the significance analysis.
#' @param shuffling_mode character, default value is "global". With this argument set of genes is selected for bootstrapping. When set to "global" all provided FC values will be used to generate random FC values for the pathway. When set to "local" FC values of pathway genes will be shuffled each time and used.
#' @param sum logical, default value is FALSE. When set to true pathway activity will be calculated via addition, when set to false then activity will be calculated via multiplication.
#' @param split logical, if true then the incoming signal will be proportionally spitted among the edges.
#' @param tmm_mode logical, default value is FALSE. When set to true specific PSF configuration will be used for the pathway activity calculation described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full
#' @param tmm_update_mode logical, default value is FALSE. When set to true specific PSF configuration will be used for the pathway activity calculation described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full
#' @param ncores integer, the number of CPUs to use for PSF calculation. Default value is 1.
#' @return Kegg collection list with evaluated PSF formulas, calculated PSF values and significance values. Each object of the list corresponds to one pathway where "psf_activities", "p_val_mat", "exp_fc" matrices are added (each column of the matrix is one sample). This list object can be further used for pathway visualization and PSF report generation.
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
  
  
  ### running psf on pathway
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
    
    ### calculating significance (need to be improved, very slow with this setup)
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



### complete implementation of specific influenced node(s) calculation (issues, single influence node, many samples with different influence node order, etc.)
#' Performs partial influence analysis which evaluates effect of each pathway node on specific node(s) of the pathway.
#' @param pathway Pathway object with calculated PSF activities with run_psf function. 
#' @param influenced_node single id of node or vector of node ids based on which partial influence will be calculated.
#' @param influence_direction in which direction node affects the signal of target node(s). Possible values c("+", "-", "any"). Default value is "any".
#' @param sample_id sample id or vector of sample ids for which partial influence will be calculated. Default value is "all".
#' @param node_combinations number of node combinations for partial influence analysis. Note: number of PSF calculations increases by the exponent of combinations(N nodes^combinations). Default value is 1.
#' @param get_influence_matrix When set to true influence matrix will be returned instead of influential nodes. Each column in influence matrix represent the log ratio of default and neutralized node(s) psf profile.
#' @param ncores integer, the number of CPUs to use for calculation. Default value is 1.
#' @return Returns list with 2 objects. "ordered_influence_nodes" matrix with ordered nodes based on their influence where each column correspond to one sample. "mean_influence_params" dataframe with mean influence parameters for each node (NULL when analysis performed for 1 sample).
#' @import parallel
#' @import graph
#' @export
run_pi <- function(pathway, influenced_node, influence_direction = "any", sample_id = "all", node_combinations = 1, get_influence_matrix = FALSE, ncores = 1) {
  
  if(!("eval.exprs" %in% names(pathway))) {
    stop("Please provide pathway with evaluated activity")
  }
  
  if(!all(influenced_node %in% rownames(pathway$psf_activities))) {
    stop("Specified id for influenced node(s) does not exist in the pathway")
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
  
  if(sample_id != "all") {
    if(!all(sample_id %in% colnames(pathway$exp_fc))) {
      stop("Please specify correct sample id(s)")
    } else {
      pathway$exp_fc <- pathway$exp_fc[,sample_id, drop = F]
      pathway$psf_activities <- pathway$psf_activities[,sample_id, drop = F]
    }
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
    } else {
      mean_influence_params <- NULL
    }
    
    ordered_influence_nodes <- Reduce(cbind, psf_influence_list)
    
    return(list(ordered_influence_nodes = ordered_influence_nodes, mean_influence_params = mean_influence_params))
  }
  
}


#' Plots the pathway with colored nodes and labels
#' @param pathway pathway object.
#' @param plot_type network visualization type. Possible values c("kegg", "visnet"). When "kegg" option is used the pathway will be plotted over kegg png image. With "visnet" option, function will plot interactive network.
#' @param plot_layout layout type of the pathway. Only works in visnet plot_type. Default value is NULL. Available values c("layout_nicely").
#' @param color_nodes type of node values to be visualized. When value type is specified pathway nodes will be color coded with expression FC values or PSF values and color legend will be added to the pathway plot. Possible values are c(NULL, "psf_activities", "exp_fc"). Default value is NULL.
#' @param multi_color_nodes when set to true nodes of the pathway will be colored with multiple colors based on the provided sample vector or groups (maximum number is 4). Available only for the kegg visualization.
#' @param sample_id name or character vector of the samples which will be used to visualize psf or exp fc values on pathway. To get averaged values across all samples, set value to "mean". Default value is "mean". 
#' @param col_data sample information for calculated PSF activities. A data.frame with two columns, first column are sample ids (colnames of PSF activity matrix), second column group names. Default value is NULL. When data.frame with group information is provided plot_pathway will generate boxplots of sink values for each group in the righ side of the pathway.
#' @param log_norm log transform PSF and expression values before color mapping. Default value is TRUE.
#' @param highlight_nodes single value of node id or a vector of ids to be highlighted in plotted pathway. Default value is NULL.
#' @param highlight_color Highlighted nodes color(s). Default values is "red".
#' @param plot_sink_values If set to TRUE then pathway activity values of terminal (sink) nodes will be plotted on the right side of the pathway (boxplots for multiple samples, barplot for single sample). Default value is FALSE. 
#' @param y_adj_text Numeric value to adjust Y position of node labels. Depending on R version this value should be adjusted to have accurate text alignment. Default value is 0.
#' @param y_adj_sink Numeric value to adjust Y position of sink node sign labels. Depending on R version this value should be adjusted to have accurate alignment. Default value is 0.
#' @param use_old_images use old kegg images(for use with curated pathway collection). Default value is FALSE.
#' @import graph
#' @import ggplot2
#' @import visNetwork
#' @import RCurl
#' @export
plot_pathway <- function(pathway, plot_type = "visnet", plot_layout = NULL,
                         color_nodes = NULL, multi_color_nodes = FALSE, sample_id = "mean", col_data = NULL, log_norm = T, 
                         highlight_nodes = NULL, highlight_color = "red", plot_sink_values = F,
                         y_adj_text = 0, y_adj_sink = 0, use_old_images = F) {
  
  ### checking highlight nodes and their border colors
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
    
    if(multi_color_nodes) {
      
      if(length(sample_id) > 1 && !is.null(col_data)) {
        stop("Multiple samples can not be plotted when coldata is provided.")
      }
      
      if(length(unique(col_data[,2])) > 4 || length(sample_id) > 4) {
        stop("Multicolor node mode is available for up to 4 groups/samples")
      }
      
      gene_nodes <- graph::nodes(pathway$graph)[which(unlist(graph::nodeData(pathway$graph, attr = "type")) == "gene")]
      
      if(length(sample_id) > 1) {
        pathway_exp_values <- Reduce(c,
                                     lapply(sample_id, function(x) {
                                       setNames(nm = paste(x, gene_nodes, sep = ";"), object = pathway$exp_fc[gene_nodes,x])
                                     }))
        
        pathway_psf_values <- Reduce(c, 
                                     lapply(sample_id, function(x) {
                                       setNames(nm = paste(x, gene_nodes, sep = ";"), pathway$psf_activities[gene_nodes,x])
                                     }))
      }
      
      if(!is.null(col_data)) {
        pathway_exp_values <- Reduce(c,
                                     lapply(unique(col_data[,2]), function(x) {
                                       setNames(nm = paste(x, gene_nodes, sep = ";"), 
                                                object = rowMeans(pathway$exp_fc[gene_nodes,col_data[which(col_data[,2] == x),1]])
                                       )
                                     }))
        
        pathway_psf_values <- Reduce(c, 
                                     lapply(unique(col_data[,2]), function(x) {
                                       setNames(nm = paste(x, gene_nodes, sep = ";"), 
                                                rowMeans(pathway$psf_activities[gene_nodes,col_data[which(col_data[,2] == x),1]])
                                       )
                                     }))   
      }
      
      
    } else {
      if(length(sample_id) == 1) {
        if(sample_id == "mean") {
          pathway_exp_values <- rowMeans(pathway$exp_fc)
          pathway_psf_values <- rowMeans(pathway$psf_activities)
        } else {
          if(!(sample_id %in% colnames(pathway$psf_activities))) {
            stop("Provided sample Ids must match with sample names of pathway activity results")
          }
          pathway_exp_values <- rowMeans(pathway$exp_fc[,sample_id,drop=F])
          pathway_psf_values <- rowMeans(pathway$psf_activities[,sample_id,drop=F])
        }
        
      } else {
        if(!all(sample_id %in% colnames(pathway$psf_activities))) {
          stop("Provided sample Ids must match with sample names of pathway activity results")
        }
        pathway_exp_values <- rowMeans(pathway$exp_fc[,sample_id,drop=F])
        pathway_psf_values <- rowMeans(pathway$psf_activities[,sample_id,drop=F])
      }
    }
    
    
    if(log_norm) {
      if(any(drop(pathway_psf_values) < 0)) {
        print("Negative values are detected skipping log transformation")
      } else {
        pathway_exp_values <- round(log(drop(pathway_exp_values)[order(drop(pathway_exp_values))] + 0.00001), digits = 5)
        pathway_psf_values <- round(log(drop(pathway_psf_values)[order(drop(pathway_psf_values))] + 0.00001), digits = 5)
      }
    } else {
      if(any(drop(pathway_psf_values) < 0)) {
        log_norm <- TRUE
      }
      pathway_exp_values <- round(drop(pathway_exp_values)[order(drop(pathway_exp_values))], digits = 5)
      pathway_psf_values <- round(drop(pathway_psf_values)[order(drop(pathway_psf_values))], digits = 5)
    }
    
  } else {
    mapping_data = FALSE
  }
  
  ### adding group coloring. Still experimental
  # if(!is.null(col_data)) {
  #   color_nodes = NULL
  # }
  
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
      
      if(multi_color_nodes) {
        node_colors$group <- sapply(node_colors$node_id, function(x) {unlist(strsplit(x, split = ";"))[1]})
        node_colors$node_id <- sapply(node_colors$node_id, function(x) {unlist(strsplit(x, split = ";"))[2]})
      } else {
        node_colors <- node_colors[which(node_colors$node_id %in% graph::nodes(pathway$graph)[which(unlist(graph::nodeData(pathway$graph, attr = "type")) != "map")]),]
        rownames(node_colors) <- node_colors$node_id
      }
      
      
    } else {
      stop("Please provide pathway with evalueated activity")
    }  
    
    if(color_nodes == "psf_activities") {
      col_legend_title = ifelse(log_norm, "Log PSF value", "PSF value")
    } else {
      col_legend_title = ifelse(log_norm, "Log FC value", "FC value")
    }
    
  } else {
    node_colors <- NULL
  }
  
  pathway_tables <- graphnel_to_df(pathway, extended = TRUE)
  
  node_graphics <- pathway_tables$node_table
  
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
    
    ### node coloring
    if(any(node_graphics$node_id %in% node_colors$node_id)) {
      if(multi_color_nodes) {
        
        group_size <- length(unique(node_colors$group))
        
        lapply(0:(length(unique(node_colors$group)) - 1), function(x) {
          
          group_node_colors <- node_colors[which(node_colors$group == unique(node_colors$group)[x+1]),]
          coloring_set <- node_graphics[unique(group_node_colors$node_id),]
          coloring_set$width <- (coloring_set$x_end - coloring_set$x_start)/group_size
          
          coloring_set$x_start <- coloring_set$x_start + coloring_set$width*x
          coloring_set$x_end <- coloring_set$x_start + coloring_set$width
          
          graphics::rect(coloring_set$x_start,
                         coloring_set$y_start - 1,
                         coloring_set$x_end,
                         coloring_set$y_end,
                         # border = coloring_set$border_color,
                         border = NA,
                         # lty = coloring_set$lty_type,
                         lwd=2,
                         col = grDevices::adjustcolor(group_node_colors$col, alpha.f = 1)
          )
          
          if(x == (length(unique(node_colors$group)) - 1)) {
            graphics::text(x = coloring_set$x,
                           y = coloring_set$y_start + y_adj_text,
                           labels = coloring_set$label,
                           col = "black", adj = c(0,0.2) + c(0.48, 1))
          }
          
        })
        
      } else {
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
        
        graphics::text(x = coloring_set$x,
                       y = coloring_set$y_start + y_adj_text,
                       labels = coloring_set$label,
                       col = node_colors$text_col, adj = c(0,0.2) + c(0.48, 1))
        
      }
    }
    
    ### adding sink node labels
    if(sum(node_graphics$sink) > 0) {
      graphics::text(x = node_graphics[node_graphics$sink,"x_end"] + 10,
                     y = node_graphics[node_graphics$sink,"y"] - 30 + y_adj_sink, cex = 3,
                     labels = rep("*", sum(node_graphics$sink)),
                     col = rep("#9ACD32", sum(node_graphics$sink)), adj = c(0,0.2) + c(0.48, 1))
    }
    
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
    
    
    ## color group nodes
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
    
    if(plot_sink_values) {
      if(mapping_data) {
        
        sink_names <- unlist(graph::nodeData(pathway$graph, pathway$sink.nodes, attr = "label"))
        
        sink_order_by_coord <- order(as.integer(unname(unlist(graph::nodeData(pathway$graph, pathway$sink.nodes, attr = ifelse(any(grepl("kegg", names(graph::nodeDataDefaults(pathway$graph)))), "kegg.gr.y", "y"))))), decreasing = T)
        
        repeating_sinks <- table(sink_names)[which(table(sink_names) > 1)]
        
        if(length(repeating_sinks) > 0) {
          for(i in 1:length(repeating_sinks)) {
            
            sink_names[which(sink_names == names(repeating_sinks[i]))] <- sapply(1:repeating_sinks[i], function(x) {paste0(paste(rep(" ", x - 1), collapse = ""), names(repeating_sinks[i]))})
            
          }
        }
        
        
        ### modify here for group separate boxplots
        psf_activities_vec <- as.vector(pathway$psf_activities[pathway$sink.nodes,])
        
        names(psf_activities_vec) <- rep(pathway$sink.nodes, ncol(pathway$psf_activities))
        
        if(log_norm) {
          if(!any(drop(pathway_psf_values) < 0)) {
            psf_activities_vec <- log(psf_activities_vec)
          }
        }
        
        psf_activity_order <- order(psf_activities_vec)
        
        psf_activities_vec <- psf_activities_vec[psf_activity_order]
        
        if(!is.null(col_data)) {
          
          rownames(col_data) <- col_data[,1]
          
          sample_group <- unlist(lapply(colnames(pathway$psf_activities), function(x) {rep(col_data[x,2], length(pathway$sink.nodes))}))
          
          sample_group <- sample_group[psf_activity_order]
          
        } else {
          sample_group <- ""
        }
        
        sink_signals <- data.frame(sink_name = sink_names[names(psf_activities_vec)], 
                                   signal = psf_activities_vec, 
                                   dot_color = color_code(values = psf_activities_vec, pal1 = pal1, pal2 = pal2, log_scale = log_norm),
                                   sample_group = sample_group,
                                   stringsAsFactors = F)
        
        sink_signals$sink_name <- factor(sink_signals$sink_name, levels = sink_names[sink_order_by_coord])
        
        if(ncol(pathway$exp_fc) == 1 | sample_id != "mean") {
          
          sink_plot <- ggplot(sink_signals, aes(x=sink_name, y=signal)) +
            geom_bar(color="black", lwd=0.2, stat="identity") +
            coord_flip() +
            theme_bw()
          
        } else {
          sink_plot <- ggplot(sink_signals, aes(x=sink_name, y=signal, fill=signal)) +
            geom_boxplot(color="black", lwd=0.2, outlier.shape=NA) +
            geom_point(color = sink_signals$dot_color) +
            geom_jitter(color = sink_signals$dot_color, width = 0.2) +
            guides(fill=FALSE) +
            coord_flip() +
            facet_wrap(~sample_group, scale="fixed", ncol = length(unique(sink_signals$sample_group))) +
            theme_bw()
        }
        
        ggsave("sink_plot.png", plot = sink_plot, device = "png", path = NULL,
               scale = 1, width = 400*length(unique(sink_signals$sample_group)), height = magick::image_info(img)$height, units = "px",
               dpi = 96, limitsize = TRUE)
        
        sink_plot_img <- image_read('sink_plot.png')
        
        file.remove('sink_plot.png')
        
        combined_img <- image_append(c(img, sink_plot_img))
        
        return(combined_img)
        
      } else {
        stop("Please provide pathway with evalueated activity")
      }
      
    } else {
      return(img)
    }
    
  } else {
    
    node_table <- pathway_tables$node_table[,c("node_id", "label", "vis_width", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title", "x", "y")]
    
    colnames(node_table) <- c("id", "label", "widthConstraint.maximum", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title", "x", "y")
    
    if(!is.null(plot_layout)) {
      node_table <- node_table[,1:11]
    }
    
    edge_table <- pathway_tables$edge_table[,c("from", "to", "color", "arrows.to.enabled", "arrows.to.type", "label", "dashes")]
    
    if(!is.null(node_colors)) {
      node_table$color.background <- unname(sapply(node_table$id, function(x) {
        if(x %in% node_colors$node_id) {
          node_colors[which(node_colors$node_id == x),"col"]
        } else {
          "#BFFFBF"
        }
      }))
      
      node_table$font.color <- unname(sapply(node_table$id, function(x) {
        if(x %in% node_colors$node_id) {
          node_colors[which(node_colors$node_id == x),"text_col"]
        } else {
          "#000000"
        }
      }))
      
      if(log_norm) {
        node_table$title <- paste(node_table$title,
                                  paste("Log exp FC", pathway_exp_values[node_table$id]),
                                  paste("Log PSF", pathway_psf_values[node_table$id]),
                                  sep = "<br>")
      } else {
        node_table$title <- paste(node_table$title,
                                  paste("Exp FC", pathway_exp_values[node_table$id]),
                                  paste("PSF", pathway_psf_values[node_table$id]),
                                  sep = "<br>")
      }
      
    } else {
      
      if(mapping_data) {
        node_table$title <- paste(node_table$title,
                                  paste("Exp FC", pathway_exp_values[node_table$id]),
                                  paste("PSF", pathway_psf_values[node_table$id]),
                                  sep = "<br>")
      }
      
    }
    
    if(!is.null(highlight_nodes)) {
      
      node_table$color.border <- unname(sapply(node_table$id, function(x) {
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
    
    node_table$image <- "unselected"
    
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
        id = as.character(max(as.integer(node_table$id)) + 1),
        label = "Color legend", 
        widthConstraint.maximum = NA, 
        font.color = "#000000", 
        size = 30, 
        font.size = 32, 
        borderWidth = 0, 
        color.border = "",
        color.background = "",
        shape = "image", 
        title = "",
        image = legend_path
      )
      
      if(is.null(plot_layout)) {
        legend_data_frame$x <- max(node_table$x)
        legend_data_frame$y <- min(node_table$y)
      }
      
      node_table <- rbind(node_table, legend_data_frame)
    }
    
    saving_pathway_name <- gsub("(", "", pathway$attrs$title, fixed = T)
    saving_pathway_name <- gsub(")", "", saving_pathway_name, fixed = T)
    saving_pathway_name <- gsub("-", "_", saving_pathway_name)
    saving_pathway_name <- gsub(" ", "_", saving_pathway_name)
    saving_pathway_name <- gsub("___", "_", saving_pathway_name)
    
    visNetwork(nodes = node_table, edges = edge_table, width = "100%", height = "800px") %>% 
      visIgraphLayout(layout = "layout_nicely") %>% 
      visExport(name = saving_pathway_name) %>%
      visInteraction(navigationButtons = TRUE, multiselect = T)
    
  }
  
}


#' Plots the TMM pathway with colored nodes and labels with interactive network.
#' @param pathway pathway object.
#' @param no_color_mode whetherto colorcode nodes based on node values or not. Default value is FALSE.
#' @param mapping_data_type type of node values to be visualized. When value type is specified pathway nodes will be color coded with expression FC values or PSF values and color legend will be added to the pathway plot. Possible values are c(NULL, "signal", "exp_fc"). Default value is signal.
#' @param log_norm log transform PSF and expression values before color mapping. Default value is TRUE.
#' @param layout layout type of the pathway. Default value is layout_nicely. Available values c("layout_nicely").
#' @import graph
#' @import visNetwork
#' @import RCurl
#' @export
plot_tmm_pathway <- function(pathway, no_color_mode = F, mapping_data_type = "signal", log_norm = TRUE, layout = "layout_nicely") {
  
  exp_values_all <- unlist(graph::nodeData(pathway$graph, attr = "expression"))[which(unlist(graph::nodeData(pathway$graph, attr = "type")) == "gene")]
  
  if(log_norm) {
    mean_exp_values <- round(log(drop(exp_values_all)[order(drop(exp_values_all))] + 0.00001), digits = 5)
  } else {
    mean_exp_values <- round(drop(exp_values_all)[order(drop(exp_values_all))], digits = 5)
  }
  
  if(no_color_mode) {
    exp_colors <- NULL
  } else {
    
    exp_colors <- color_code(values = mean_exp_values, pal1 = pal1, pal2 = pal2, log_scale = log_norm)
    
    exp_colors = data.frame(node_id = names(mean_exp_values)[c(which(mean_exp_values <= 0), which(mean_exp_values > 0))], 
                            col = exp_colors,
                            exps = mean_exp_values[c(which(mean_exp_values <= 0), which(mean_exp_values > 0))],
                            text_col = unname(sapply(exp_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                            stringsAsFactors = F
    )
  }
  
  
  signal_values_all <- unlist(graph::nodeData(pathway$graph, attr = "signal"))[which(unlist(graph::nodeData(pathway$graph, attr = "type")) != "map")]
  
  if(log_norm) {
    mean_signal_values <- round(log(signal_values_all[order(signal_values_all)] + 0.00001), digits = 5)
  } else {
    mean_signal_values <- round(signal_values_all[order(signal_values_all)], digits = 5)
  }
  
  if(no_color_mode) {
    psf_colors <- NULL
  } else {
    psf_colors <- color_code(values = mean_signal_values, pal1 = pal1, pal2 = pal2, log_scale = log_norm)
    
    psf_colors = data.frame(node_id = names(mean_signal_values)[c(which(mean_signal_values <= 0), which(mean_signal_values > 0))], 
                            col = psf_colors,
                            signals = mean_signal_values[c(which(mean_signal_values <= 0), which(mean_signal_values > 0))],
                            text_col = unname(sapply(psf_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                            stringsAsFactors = F
    )
  }
  
  if(mapping_data_type == "signal") {
    color.genes <- psf_colors
    color_bar_lims <- range(mean_signal_values)
    col_legend_title = ifelse(log_norm, "Log PSF value", "PSF value")
  } else {
    color.genes <- exp_colors
    color_bar_lims <- range(mean_exp_values)
    col_legend_title = ifelse(log_norm, "Log FC value", "FC value")
  }
  
  
  ### network graphics building
  
  entrez_id <- unname(sapply(pathway$graph@nodes, function(y) {
    ifelse(unlist(graph::nodeData(pathway$graph, y, attr = "genes")) == 0,
           as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))), 
           paste(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), collapse = ","))
  }))
  
  hover_name <- unname(sapply(pathway$graph@nodes, function(y) {
    ifelse(unlist(graph::nodeData(pathway$graph, y, attr = "genes")) == 0,
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
  
  node_coords <- data.frame(
    node_type = as.character(unlist(graph::nodeData(pathway$graph, attr = "type"))),
    psf_function = as.character(unlist(graph::nodeData(pathway$graph, attr = "psf_function"))),
    gr_name = as.character(unlist(graph::nodeData(pathway$graph, attr = "label"))),
    node_id = names(graph::nodeData(pathway$graph, attr = "type")),
    entrez_id = entrez_id, 
    hover_name = paste(hover_name, as.character(unlist(graph::nodeData(pathway$graph, attr = "type"))), sep = "_"),
    sink = (names(graph::nodeData(pathway$graph, attr = "type")) %in% pathway$sink.nodes), 
    x = unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.x")),
    y = unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.y")),
    stringsAsFactors = F 
  )
  
  kegg_arrows_type <- c("simple", "simple", "simple", "simple", "T", "simple", "simple", "simple", "T", "simple", "", "simple", "simple", "T", "simple", "simple")
  names(kegg_arrows_type) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
  
  line_col <- c("red", "red", "red", "red", "blue","red", "red", "red", "blue", "red", "red", "red", "red", "blue", "red", "red")
  names(line_col) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
  
  
  if(length(names(pathway$graph@edgeData@data)) > 0) {
    splitted_interactions <- strsplit(names(pathway$graph@edgeData@data), split = "|", fixed = T)
    from <- as.character(sapply(splitted_interactions, "[[", 1))
    to <- as.character(sapply(splitted_interactions, "[[", 2))
  } else {
    edge_coords < NULL
  }
  
  from_index <- match(from, node_coords$node_id)
  to_index <- match(to, node_coords$node_id)
  
  col <- line_col[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype1")))]
  
  edge_coords <- data.frame(
    from = from, to = to,
    arr.type = kegg_arrows_type[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "subtype1")))],
    col = col, stringsAsFactors = F
  )
  
  if(!is.null(color.genes)) {
    color <- unname(sapply(node_coords$node_id, function(x) {
      if(x %in% color.genes$node_id) {
        color.genes[which(color.genes$node_id == x),"col"]
      } else {
        "#BFFFBF"
      }
    }))
    
    if(log_norm) {
      hover_string <- paste(node_coords$hover_name, 
                            paste("Node id", node_coords$node_id),
                            paste("PSF function", node_coords$psf_function),
                            paste("Log exp FC", unname(sapply(node_coords$node_id, function(x) {exp_colors[which(exp_colors$node_id == x),"exps"]}))), 
                            paste("Log PSF", unname(sapply(node_coords$node_id, function(x) {psf_colors[which(psf_colors$node_id == x),"signals"]}))), 
                            sep = "<br>")
    } else {
      hover_string <- paste(node_coords$hover_name, 
                            paste("Node id", node_coords$node_id),
                            paste("PSF function", node_coords$psf_function),
                            paste("Exp FC", unname(sapply(node_coords$node_id, function(x) {exp_colors[which(exp_colors$node_id == x),"exps"]}))), 
                            paste("PSF", unname(sapply(node_coords$node_id, function(x) {psf_colors[which(psf_colors$node_id == x),"signals"]}))), 
                            sep = "<br>")
    }
    
    
    
    font_color <- unname(sapply(node_coords$node_id, function(x) {
      if(x %in% color.genes$node_id) {
        color.genes[which(color.genes$node_id == x),"text_col"]
      } else {
        "#000000"
      }
    }))
    
  } else {
    
    hover_string <- paste(node_coords$hover_name, 
                          paste("Node id", node_coords$node_id),
                          paste("PSF function", node_coords$psf_function),
                          paste("Exp FC", exp_values_all[node_coords$node_id]),
                          paste("PSF", signal_values_all[node_coords$node_id]),
                          sep = "<br>")
    
    color_converter <- setNames(object = c("#a7d1dfff", "#a7d1dfff", "#e3d161ff", "#e3d161ff", "#e3d161ff", "#9a9a9aff", "#9a9a9aff", "#a7d1dfff"),
                                nm = c("gene",  "C", "BP", "BP.0", "BP.t", "L-max", "L-min", "S"))
    
    color <- color_converter[node_coords$node_type]
    
    font_color <- rep("#000000", nrow(node_coords))
    
  }
  
  # border_color_converter <- setNames(object = c("#ffff00", "#66ff33", "#BFFFBF"), nm = c("human", "SARS-COV-2", "---"))
  
  border_color <- unname(sapply(1:length(node_coords$sink), function(x) {
    if(node_coords$sink[x]) {
      "#0099cc"
    } else {
      color[x]
    }
  }))
  
  node_type_and_shape <- setNames(object = c("box", "ellipse", "ellipse", "ellipse", "ellipse", "triangle", "triangle", "ellipse"),
                                  nm = c("gene",  "C", "BP", "BP.0", "BP.t", "L-max", "L-min", "S"))
  
  
  nodes <- data.frame(id = node_coords$node_id,
                      image = rep("unselected", nrow(node_coords)),
                      label = node_coords$gr_name,
                      shape = node_type_and_shape[node_coords$node_type],
                      color.background = color, 
                      color.border = border_color, 
                      # color.highlight = color,
                      borderWidth = 2,
                      title = hover_string,
                      font.size = rep(22, nrow(node_coords)), size = 25,
                      font.color = font_color, row.names = node_coords$node_id,
                      x = node_coords$x,
                      y = node_coords$y
  )
  
  if(!is.null(color.genes)) {
    
    # legend_img <- magick::image_device(width = 381, height = 280)
    legend_img <- magick::image_device(width = 480, height = 480)
    plot.new()
    
    color_legend_maker(x = 0.05, y = 0, leg = 0.9, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                       lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    
    
    magick::image_write(magick::image_trim(legend_img, fuzz = 0), path = "~/legend.png")
    # magick::image_write(magick::image_crop(legend_img, geometry = "400x30+45+370"), path = "~/legend.png")
    # magick::image_write(magick::image_border(magick::image_crop(legend_img, geometry = "420x30+45+370", repage = F), color = "blue"), path = "~/legend.png")
    
    legend_path <- paste('data:image/png;base64', RCurl::base64Encode(readBin('~/legend.png', 'raw', file.info('~/legend.png')[1, 'size']), 'txt'), sep = ',')
    
    legend_data_frame <- data.frame(
      id = as.character(max(as.integer(node_coords$node_id)) + 1),
      image = legend_path,
      label = "Color legend", shape = "image", color.background = "",
      color.border = "", borderWidth = 0, title = "",
      font.size = 32, size = 25, font.color = "#000000",
      x = max(node_coords$x),
      y = min(node_coords$y) - 10
    )
    
    nodes <- rbind(nodes, legend_data_frame)
  }
  
  arrows_type <- c("arrow","bar")
  names(arrows_type) <- c("simple", "T")
  
  edges <- data.frame(from = edge_coords$from, to = edge_coords$to,
                      color = edge_coords$col,
                      arrows.to.enabled = rep(TRUE, length(edge_coords$col)),
                      arrows.to.type = arrows_type[edge_coords$arr.type]
  )
  
  lnodes = data.frame(label = c("gene", "Isoform", "Complex/state", "Process", "Linker(min/max/sum)"),
                      shape = c("box", "box", "ellipse", "ellipse", "triangle"),
                      sizes = 25,
                      color.background = c("#a7d1dfff", "#a7d1dfff", "#a7d1dfff", "#e3d161ff", "#9a9a9aff"),
                      color.border = ""
  )
  
  ledges = data.frame(label = c("Activation", "Inhibition"), 
                      arrows.to.enabled = TRUE, 
                      arrows.to.type = c("arrow", "bar"), 
                      color = c("red", "blue"))
  
  network <- visNetwork(nodes = nodes, edges = edges, width = "100%", height = "800px") %>% 
    visIgraphLayout_new(layout = layout) %>% 
    visInteraction(navigationButtons = TRUE, multiselect = T) %>%
    visLegend(addEdges = ledges, addNodes = lnodes, useGroups = FALSE, width = 0.1, stepY = 60, position = "right", zoom = FALSE)
  
  return(network)
  
}




#' Generates pdf report with colored pathways and plots
#' @param psf_list Output from run_psf function
#' @param folder_name name of the folder where pdf report(s) will be generated
#' @param plot_type Pathway visualization type. Possible values are c("kegg", "visnet"). Default value is "kegg".
#' @param log_norm log transform PSF values before color mapping. Default value is TRUE.
#' @param y_adj_text Numeric value to adjust Y position of node labels. Depending on R version this value should be adjusted to have accurate text alignment. Default value is 0.
#' @param y_adj_sink Numeric value to adjust Y position of sink node sign labels. Depending on R version this value should be adjusted to have accurate alignment. Default value is 0.
#' @param use_old_images use olde kegg images(for use with curated pathway collection)
#' @import gplots
#' @import ggplot2
#' @import ggrepel
#' @import magick
#' @import grDevices
#' @export
generate_psf_report <- function(psf_list, folder_name, plot_type = "kegg", log_norm = TRUE, y_adj_text = 0, y_adj_sink = 0, use_old_images = F) {
  
  dir.create(folder_name)
  
  lapply(psf_list, function(x) {
    
    magick::autoviewer_disable()
    
    print(paste0("Generating report for ", x$attrs$title, " pathway"))
    
    pathway_plot <- plot_pathway(pathway = x, plot_type = plot_type, color_nodes = "psf_activities", 
                                 sample_id = "mean", log_norm = log_norm, plot_sink_values = T, 
                                 y_adj_text = y_adj_text, y_adj_sink = y_adj_sink,
                                 use_old_images = use_old_images)
    
    plots <- magick::image_graph(width = 1000, height = 800, res = 96)
    
    if(ncol(x$exp_fc) > 1) {
      if(log_norm) {
        sink_values <- log(x$psf_activities[x$sink.nodes,])
      } else {
        sink_values <- x$psf_activities[x$sink.nodes,]
      }
      gplots::heatmap.2(sink_values, col = c(pal1(10), pal2(10)), trace = "none",
                        Rowv = FALSE, Colv = FALSE, ylab = "Sink values", main = paste0(x$attrs$title, " Sink PSF values"),
                        margin = c(10, 8), keysize = 1, key.title = "PSF log value",
                        dendrogram = 'none')
    }
    dev.off()
    
    if("p_val_mat" %in% names(x)) {
      
      ### generating volcano plot
      sink_names <- unlist(graph::nodeData(x$graph, x$sink.nodes, attr = "label"))
      
      psf_activities_vec <- as.vector(x$psf_activities[x$sink.nodes,])
      
      volcano_table <- data.frame(sink_names = unlist(lapply(colnames(x$psf_activities), function(x) {paste(sink_names, x, sep = "_")})), 
                                 log_psf = log(psf_activities_vec), 
                                 pvalue = as.vector(x$p_val_mat[x$sink.nodes,]) + 0.00001,
                                 stringsAsFactors = F)
      
      ## detecting up and down regulated significan sink psf values
      volcano_table$diffpsf <- "NO"
      volcano_table$diffpsf[volcano_table$log_psf > 0.6 & volcano_table$pvalue < 0.05] <- "UP"
      volcano_table$diffpsf[volcano_table$log_psf < -0.6 & volcano_table$pvalue < 0.05] <- "DOWN"
      
      ## adding labels for signifcant sinks
      volcano_table$psflabel <- ""
      volcano_table$psflabel[volcano_table$diffpsf != "NO"] <- volcano_table$sink_names[volcano_table$diffpsf != "NO"]
      
      
      volcano_plot <- ggplot(data=volcano_table, aes(x=log_psf, y=-log10(pvalue), col=diffpsf, label=psflabel)) +
        geom_point() + 
        theme_minimal() +
        ggtitle("PSF with p values") +
        theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) +
        geom_text_repel(max.overlaps = Inf) +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-0.6, 0.6), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")
      
      ggsave("volcano_plot.png", plot = volcano_plot, device = "png", path = NULL,
             scale = 1, width = 1000, height = 800, units = "px",
             dpi = 96, limitsize = TRUE)
      
      volcano_plot_img <- image_read('volcano_plot.png')
      
      file.remove('volcano_plot.png')
      
      plots <- image_append(c(plots, volcano_plot_img))
      
    }
    
    magick::image_write(c(pathway_plot, plots), format = "pdf", path = paste0(folder_name, "/", gsub("[- ]", "_", x$attrs$title), ".pdf"), quality = 300)
    
  })
  
}

#' Converts graphNEL object to 2 data frames(node_table, edge_table)
#' @param pathway pathway list object generated with package.
#' @param extended set to TRUE to get more detailed data frames with graphical information. Default value is FALSE.
#' @import graph
#' @export
graphnel_to_df <- function(pathway, extended = FALSE) {
  
  if(length(pathway$graph@edgeData@data) > 0) {
    splitted_interactions <- strsplit(names(pathway$graph@edgeData@data), split = "|", fixed = T)
    from <- as.character(sapply(splitted_interactions, "[[", 1))
    to <- as.character(sapply(splitted_interactions, "[[", 2))
  }
  
  if(any(grepl("kegg", names(graph::nodeDataDefaults(pathway$graph))))) {
    attribute_converter <- setNames(nm = c("kegg.id", "kegg.gr.x", "kegg.gr.y", "kegg.gr.width", "kegg.gr.height"), object = c("kegg.id", "kegg.gr.x", "kegg.gr.y", "kegg.gr.width", "kegg.gr.height"))
  } else {
    attribute_converter <- setNames(nm = c("kegg.id", "kegg.gr.x", "kegg.gr.y", "kegg.gr.width", "kegg.gr.height"), object = c("node_id", "x", "y", "node_width", "node_height"))
  }
  
  # if(pathway_source == "editor") {
  #   attribute_converter <- setNames(nm = c("kegg.id", "kegg.gr.x", "kegg.gr.y", "kegg.gr.width", "kegg.gr.height"), object = c("node_id", "x", "y", "node_width", "node_height"))
  # } else {
  #   attribute_converter <- setNames(nm = c("kegg.id", "kegg.gr.x", "kegg.gr.y", "kegg.gr.width", "kegg.gr.height"), object = c("kegg.id", "kegg.gr.x", "kegg.gr.y", "kegg.gr.width", "kegg.gr.height"))
  # }
  
  
  if(extended) {
    hover_name <- unname(sapply(pathway$graph@nodes, function(y) {
      ifelse(length(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))) < 1,
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
    
    hover_name <- paste(hover_name,
                        paste("Node id", as.character(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.id"])))), sep = "<br>")
    
    
    node_shapes <- unname(sapply(as.character(unlist(graph::nodeData(pathway$graph, attr = "type"))), function(x) {
      if(grepl("compound",x)) {
        "dot"
      } else {
        "box"
      }
    }))
    
    
    color <- unname(sapply((unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.id"])) %in% pathway$sink.nodes), function(x) {
      if(x) {
        "#0099cc"
      } else {
        "#BFFFBF"
      }
    }))
    
    size <- unname(sapply(as.character(unlist(graph::nodeData(pathway$graph, attr = "type"))), function(x) {
      if(grepl("compound",x)) {
        10
      } else {
        25
      }
    }))
    
    vis_width <- unname(sapply(pathway$graph@nodes, function(x) {
      if(grepl("map", unlist(graph::nodeData(pathway$graph, x, attr = "type")))) {
        as.numeric(unlist(graph::nodeData(pathway$graph, x, attr = attribute_converter["kegg.gr.width"])))*2
      } else {
        NA
      }
    }))
    
    
    node_table <- data.frame(node_id = as.character(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.id"]))),
                             label = as.character(unlist(graph::nodeData(pathway$graph, attr = "label"))),
                             component_id_s = unname(sapply(pathway$graph@nodes, function(y) {
                               paste(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), collapse = ",")
                             })),
                             type = as.character(unlist(graph::nodeData(pathway$graph, attr = "type"))),
                             psf_function = as.character(unlist(graph::nodeData(pathway$graph, attr = "psf_function"))),
                             expression = as.numeric(unlist(graph::nodeData(pathway$graph, attr = "expression"))),
                             signal = as.numeric(unlist(graph::nodeData(pathway$graph, attr = "signal"))),
                             sink = (unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.id"])) %in% pathway$sink.nodes),
                             x_start = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.x"]))) - as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.width"])))*0.5, 
                             y_start = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.y"]))) - as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.height"])))*0.5,
                             x_end  = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.x"]))) + as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.width"])))*0.5,
                             y_end =  as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.y"]))) + as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.height"])))*0.5,
                             node_width = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.width"]))),
                             node_height = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.height"]))),
                             title = hover_name,
                             shape = node_shapes,
                             color.background = color,
                             color.border = "#BFFFBF",
                             borderWidth = 2,
                             font.size = 22,
                             size = size,
                             font.color = "#000000",
                             vis_width = vis_width,
                             x = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.x"]))),
                             y = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.y"]))),
                             existence = "exist",
                             change_info = "no_change",
                             data_source = "kegg",
                             row.names = as.character(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.id"]))),
                             stringsAsFactors = F
    )
    
    kegg_arrows_type <- c("arrow", "arrow", "arrow", "arrow", "bar", "arrow", "arrow", "arrow", "bar", "arrow", "", "arrow", "arrow", "bar", "arrow", "arrow")
    names(kegg_arrows_type) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
    
    line_col <- c("red", "red", "red", "red", "blue","red", "red", "red", "blue", "red", "red", "red", "red", "blue", "red", "red")
    names(line_col) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
    
    if(length(pathway$graph@edgeData@data) > 0) {
      edge_table <- data.frame(id = paste0(from, "|", to),
                               from = from,
                               to = to,
                               color = line_col[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "subtype1")))],
                               arrows.to.enabled = TRUE,
                               arrows.to.type = kegg_arrows_type[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "subtype1")))],
                               label = "",
                               dashes = unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "subtype2"))) == "indirect effect" | unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "subtype1"))) == "indirect effect",
                               type = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "type")),
                               subtype1 = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype1")),
                               subtype2 = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype2")),
                               state = "",
                               weight = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "weight")),
                               existence = "exist",
                               change_info = "no_change",
                               data_source = "kegg",
                               stringsAsFactors = F
      )
    } else {
      edge_table <- NULL
    }
    
    
  } else {
    
    node_table <- data.frame(node_id = as.character(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.id"]))),
                             label = as.character(unlist(graph::nodeData(pathway$graph, attr = "label"))),
                             component_id_s = unname(sapply(pathway$graph@nodes, function(y) {
                               ifelse(is.null(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),
                                      as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))), 
                                      paste(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), collapse = ","))
                             })),
                             type = as.character(unlist(graph::nodeData(pathway$graph, attr = "type"))),
                             psf_function = as.character(unlist(graph::nodeData(pathway$graph, attr = "psf_function"))),
                             expression = as.numeric(unlist(graph::nodeData(pathway$graph, attr = "expression"))),
                             signal = as.numeric(unlist(graph::nodeData(pathway$graph, attr = "signal"))),
                             x = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.x"]))),
                             y = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.y"]))),
                             existence = "exist",
                             change_info = "no_change",
                             data_source = "kegg",
                             stringsAsFactors = F
    )
    
    if(length(pathway$graph@edgeData@data) > 0) {
      edge_table <- data.frame(from = from,
                               to = to,
                               type = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "type")),
                               subtype1 = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype1")),
                               subtype2 = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype2")),
                               state = "",
                               weight = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "weight")),
                               stringsAsFactors = F
      )
    } else {
      edge_table <- NULL  
    }
    
    
  }
  
  return(list(node_table = node_table, edge_table = edge_table))
  
}

#' Converts 2 data frames(node_table, edge_table) to graphNEL object
#' @param node_table node data frame with specified columns. Mandatory columns: node_id, label, component_id_s, type. Non mandatory columns: psf_function, x, y.
#' @param edge_table edge data frame with specified columns. Mandatory columns: from, to, type. Non mandatory columns subtype, state, weight.
#' @importFrom "igraph" "as_graphnel"
#' @importFrom "igraph" "graph_from_data_frame"
#' @import graph
#' @export
df_to_graphnel <- function(node_table, edge_table) {
  
  if(sum(c("node_id", "label", "component_id_s", "type") %in% colnames(node_table)) != 4) {
    stop("Please provide all mandatory columns of the node table.")
  }
  
  if(sum(c("from", "to", "type") %in% colnames(edge_table)) != 3) {
    stop("Please provide all mandatory columns of the edge table.")
  }
  
  # creating graphNEL object
  g <- as_graphnel(
    graph_from_data_frame(edge_table[,c("from", "to")], directed = T, vertices = node_table$node_id)
  )
  
  ### setting node default values
  graph::nodeDataDefaults(g, attr="genes") <- 0
  graph::nodeDataDefaults(g, attr="expression") <- 1
  graph::nodeDataDefaults(g, attr="signal") <- 1
  graph::nodeDataDefaults(g, attr="type") <- "gene"
  graph::nodeDataDefaults(g, attr="label") <- ""
  graph::nodeDataDefaults(g, attr="components") <- "NA"
  graph::nodeDataDefaults(g, attr="psf_function") <- "mean"
  graph::nodeDataDefaults(g, attr="x") <- NA
  graph::nodeDataDefaults(g, attr="y") <- NA
  graph::nodeDataDefaults(g, attr="node_width") <- NA
  graph::nodeDataDefaults(g, attr="node_height") <- NA
  graph::nodeDataDefaults(g, attr="node_id") <- "0"
  graph::nodeDataDefaults(g, attr = "existence") <- "exist"
  graph::nodeDataDefaults(g, attr = "change_info") <- "no_change"
  graph::nodeDataDefaults(g, attr = "data_source") <- "kegg"
  
  ### filling mandatory node attributes
  graph::nodeData(g, attr = "genes") <- strsplit(node_table$component_id_s, split = ",")
  graph::nodeData(g, attr = "label") <- node_table$label
  graph::nodeData(g, attr = "type") <- node_table$type
  graph::nodeData(g, attr = "node_id") <- node_table$node_id
  
  ### filling non mandatory node attributes
  node_additional_attrs <- setdiff(colnames(node_table), c("node_id", "label", "component_id_s", "type", "sink", "x_start", "y_start", "x_end", "y_end", "vis_width", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title"))
  
  for (i in node_additional_attrs) {
    graph::nodeData(g, attr = i) <- node_table[,i]
  }
  
  graph::edgeDataDefaults(g, attr = "impact") <- 1
  graph::edgeDataDefaults(g, attr = "weight") <- 1
  graph::edgeDataDefaults(g, attr = "type") <- NA
  graph::edgeDataDefaults(g, attr = "subtype1") <- NA
  graph::edgeDataDefaults(g, attr = "subtype2") <- NA
  graph::edgeDataDefaults(g, attr = "state") <- ""
  graph::edgeDataDefaults(g, attr = "existence") <- "exist"
  graph::edgeDataDefaults(g, attr = "change_info") <- "no_change"
  graph::edgeDataDefaults(g, attr = "data_source") <- "kegg"
  
  
  edge_additional_attrs <- setdiff(colnames(edge_table), c("id", "from", "to", "color", "arrows.to.enabled", "arrows.to.type", "label", "dashes"))
  
  for(i in edge_additional_attrs) {
    for(j in 1:nrow(edge_table)) {
      graph::edgeData(g, from = edge_table[j,"from"], to = edge_table[j,"to"], attr = i) <- edge_table[j,i]
    }
  }
  
  pathway <- list(graph = g)
  pathway$sink.nodes <- psf::determine.sink.nodes(pathway)
  pathway$order <-  psf::order.nodes(pathway$graph)
  pathway$graph <- psf::set.edge.impacts(pathway$graph)
  return(pathway)
  
}

#' Returns vector of node ids which do not have incoming edges but only outgoing.
#' @param pathway pathway list object generated with package.
#' @export
determine.input.nodes <- function(pathway) {
  dfs <- graphnel_to_df(pathway)
  return(setdiff(dfs$edge_table$from, dfs$edge_table$to))
}


#' Performs pathway activity analysis of Spatial transcriptomics data and subsequent clustering with Seurat clustering. Input data is a Seurat object.
#' @param spatial_obj Seurat object of spatial transcriptomic dataset.
#' @param pathway_collection pathway collection list which will be used of PSF analysis
#' @param gene_symbol_to_entrez gene id conversion vector where names are the gene identifiers corresponding to count matrix row names and objects are entrez gene ids.
#' @param nthreads integer, the number of CPUs to use for calculation. Default value is 1.
#' @param return_only_shiny_vars when set tot TRUE list of variables will be return which can be used to lunch PSF spatial browser shiny app. Default value is TRUE.
#' @export
spatial_psf_analysis <- function(spatial_obj, pathway_collection, gene_symbol_to_entrez, nthreads = 1, return_only_shiny_vars = TRUE) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("This function requires the 'Seurat'. Please install it using install.packages('Seurat').")
  }
  
  
  ### Exp data Normalization
  cat("Exp SCT transform (normalization)\n")
  spatial_obj <- SCTransform(spatial_obj, assay = "Spatial", verbose = FALSE)
  
  ### Exp clustering
  cat("Exp clustering\n")
  spatial_obj <- RunPCA(spatial_obj, assay = "SCT", verbose = FALSE)
  spatial_obj <- FindNeighbors(spatial_obj, reduction = "pca", dims = 1:30)
  spatial_obj <- FindClusters(spatial_obj, verbose = FALSE, resolution = 0.8)
  spatial_obj <- RunUMAP(spatial_obj, reduction = "pca", dims = 1:30)
  
  
  ### Exp fold change calculation
  spatial_sct_matrix <- as.matrix(Seurat::Assays(spatial_obj, slot = "SCT")@counts)
  spatial_sct_matrix_fc <- (spatial_sct_matrix + 1)/rowMeans(spatial_sct_matrix + 1)
  
  rownames(spatial_sct_matrix_fc) <- gene_symbol_to_entrez[rownames(spatial_sct_matrix_fc)]
  spatial_sct_matrix_fc <- spatial_sct_matrix_fc[which(!is.na(rownames(spatial_sct_matrix_fc))),]
  
  ### Pathway signal flow analysis
  cat("Pathway activity calculation\n")
  spatial_psf <- run_psf(entrez.fc = spatial_sct_matrix_fc, kegg.collection = pathway_collection, 
                         calculate.significance = F, ncores = nthreads)
  
  
  spatial_psf_mat <- Reduce(rbind, 
                            lapply(names(spatial_psf), function(x) {
                              print(x)
                              pathway <- spatial_psf[[x]]
                              
                              activity_mat <- pathway$psf_activities[pathway$sink.nodes,,drop = F]
                              
                              rownames(activity_mat) <- paste0(rownames(activity_mat), "; ", x)
                              
                              activity_mat
                            })   
  )
  
  
  spatial_psf_obj <- CreateSeuratObject(counts = spatial_psf_mat, meta.data = spatial_obj@meta.data, assay = "Spatial")
  
  cat("Patwahy activity clustering\n")
  ### Patwhay activity data normalization
  # spatial_psf_obj <- NormalizeData(spatial_psf_obj)
  
  ### Pathway activity clustering
  spatial_psf_obj <- FindVariableFeatures(spatial_psf_obj)
  spatial_psf_obj <- ScaleData(spatial_psf_obj)
  spatial_psf_obj <- RunPCA(spatial_psf_obj, verbose = FALSE, approx=TRUE)
  spatial_psf_obj <- FindNeighbors(spatial_psf_obj, dims = 1:30)
  spatial_psf_obj <- FindClusters(spatial_psf_obj, resolution = 0.8, verbose = FALSE)
  spatial_psf_obj <- RunUMAP(spatial_psf_obj, dims = 1:30)
  
  spatial_psf_obj@images <- spatial_obj@images
  
  idents <- as.character(sort(unique(as.numeric(as.character(Idents(spatial_psf_obj))))))
  
  ### Finding Cluster specific markers
  cat("Identification of cluster specific features\n")
  ident_markers <- lapply(idents, function(x) {
    markers <- FindMarkers(spatial_psf_obj, ident.1 = x, logfc.threshold = 0.25, test.use = "wilcox")
    
    if("image" %in% names(pathway_collection[[1]]$attrs)) {
      sink_annot <- kegg_sink_to_process[gsub("-", "_", rownames(markers)), c("Pathway_name", "Sink", "Process", "Pathway")]
      colnames(sink_annot)[4] <- "Dowstream_pathway"
      cbind(markers, sink_annot)
    } else {
      markers$Pathway_name <- sapply(rownames(markers), function(y) {gsub("-", "_", unlist(strsplit(y, split = "; "))[2])})
      markers$Sink <- sapply(rownames(markers), function(y) {
        graph::nodeData(pathway_collection[[gsub("-", "_", unlist(strsplit(y, split = "; "))[2])]]$graph, n = unlist(strsplit(y, split = "; "))[1], attr = "label")
      })
    }
    
  })
  
  names(ident_markers) <- idents
  
  if(return_only_shiny_vars) {
    return(list(spatial_psf_collection = spatial_psf, psf_ident_markers = ident_markers, psf_mat = spatial_psf_mat, 
                spatial_image = spatial_psf_obj@images$slice1@image, coords = GetTissueCoordinates(object = spatial_psf_obj, scale = "lowres"), 
                meta.data = spatial_psf_obj@meta.data))
  } else {
    return(list(spatial_obj = spatial_obj, spatial_psf_obj = spatial_psf_obj, spatial_psf_collection = spatial_psf,
                psf_ident_markers = ident_markers))
  }
  
}


#' Launches PSF spatial browser Shiny application.
#' @param psf_saptial_results list object generated by spatial_psf_analysis function.
#' @importFrom "magick" "image_draw"
#' @importFrom "graphics" "rect"
#' @importFrom "graphics" "text"
#' @importFrom "grDevices" "adjustcolor"
#' @importFrom "magick" "image_info"
#' @importFrom "shape" "Arrows"
#' @importFrom "grDevices" "col2rgb"
#' @importFrom "magick" "autoviewer_disable"
#' @importFrom "magick" "image_device"
#' @importFrom "magick" "image_write"
#' @importFrom "magick" "image_trim"
#' @importFrom "RCurl" "base64Encode"
#' @importFrom "plotly" "plot_ly"
#' @importFrom "plotly" "layout"
#' @importFrom "plotly" "plotlyOutput"
#' @importFrom "magick" "image_read"
#' @import visNetwork
#' @export
run_psf_spatial_browser <- function(psf_saptial_results) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("This function requires the 'Seurat'. Please install it using install.packages('Seurat').")
  }
  
  psf_subset <- function(psf_collection, pathway_name, sample_list) {
    
    pathway <- psf_collection[[pathway_name]]
    
    pathway$psf_activities <- pathway$psf_activities[,sample_list, drop = F]
    pathway$exp_fc <- pathway$exp_fc[,sample_list, drop = F]
    
    return(pathway)
  }
  
  plot_kegg_pathway <- function(graphnel_df, group_graphics, pathway_image, 
                                node_colors = NULL, custom_edge_mapping = FALSE, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
                                col_legend_title = NULL, color_bar_lims, y_adj_sink = 3, node_fill_opacity = 0, present_node_modifications = FALSE, removed_nodes = NULL, y_adj_text = 3, custom_color_scale = NULL
  ) {
    
    ### checking highlight nodes and their border colors
    if(length(highlight_color) > 1) {
      if(length(highlight_color) == length(highlight_nodes)) {
        highlight_color_vector <- setNames(object = highlight_color, nm = highlight_nodes)
      } else {
        stop("Error: highlighted nodes and their colors must be in the same length")
      }
    } else {
      highlight_color_vector <- setNames(object = rep(highlight_color, length(highlight_nodes)), nm = highlight_nodes)
    }
    
    img <- magick::image_draw(pathway_image)
    
    node_graphics <- graphnel_df$node_table
    
    ### node exp coloring
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
      
      graphics::text(x = coloring_set$x,
                     y = coloring_set$y_start + y_adj_text,
                     labels = coloring_set$label,
                     col = node_colors$text_col, adj = c(0,0.2) + c(0.48, 1))
      
    }
    
    ### adding sink node labels
    graphics::text(x = node_graphics[node_graphics$sink,"x_end"] + 10,
                   y = node_graphics[node_graphics$sink,"y"] - 30 + y_adj_sink, cex = 3,
                   labels = rep("*", sum(node_graphics$sink)),
                   col = rep("#9ACD32", sum(node_graphics$sink)), adj = c(0,0.2) + c(0.48, 1))
    
    if(!is.null(highlight_nodes)) {
      highlight_set <- node_graphics[which(node_graphics[,"node_id"] %in% highlight_nodes),]
      
      rect( highlight_set$x_start, 
            highlight_set$y_start, 
            highlight_set$x_end, 
            highlight_set$y_end, 
            border = highlight_color_vector[highlight_set$node_id], lty = "solid", lwd=2, col = adjustcolor( "#a3297a", alpha.f = node_fill_opacity))
    }
    
    
    ### scale color bar
    if(!is.null(node_colors)) {
      if(is.null(custom_color_scale)) {
        psf:::color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = c(psf:::pal1(10), psf:::pal2(10)), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                                 lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
      } else {
        psf:::color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = c(custom_color_scale[1:50], custom_color_scale[51:100]), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                                 lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
      }
    }
    
    text(x = c(magick::image_info(img)$width - 88, magick::image_info(img)$width - 30),
         y = c(70, 65 + y_adj_sink - 6), cex = c(1.5, 3),
         labels = c("Sink node", "*"),
         col = c("#000000", "#9ACD32"), adj = c(0,0.2) + c(0.48, 1))
    
    
    ## color grop nodes
    if(length(group_graphics) > 0 ) {
      lapply(group_graphics, function(z) {
        rect( z$kegg.gr.x-z$kegg.gr.width*0.5, 
              z$kegg.gr.y+z$kegg.gr.height*0.5, 
              z$kegg.gr.x+z$kegg.gr.width*0.5, 
              z$kegg.gr.y-z$kegg.gr.height*0.5, 
              border = "yellow", lty = "dashed", lwd=2)
      })
    }
    
    if(present_node_modifications) {
      added_nodes <- node_graphics[which(node_graphics$existence == "added"),]
      if(nrow(added_nodes) > 0) {
        graphics::rect(added_nodes$x_start,
                       added_nodes$y_start - 1,
                       added_nodes$x_end,
                       added_nodes$y_end,
                       # border = coloring_set$border_color,
                       border = rep("black", nrow(added_nodes)),
                       # lty = coloring_set$lty_type,
                       lwd=1,
                       col = rep("#BFFFBF", nrow(added_nodes))
        )
        
        graphics::text(x = added_nodes$x,
                       y = added_nodes$y_start + y_adj_text,
                       labels = added_nodes$label,
                       col = rep("black", nrow(added_nodes)), adj = c(0,0.2) + c(0.48, 1))
      }
      
      if(!is.null(removed_nodes)) {
        graphics::rect(removed_nodes$x_start,
                       removed_nodes$y_start - 1,
                       removed_nodes$x_end,
                       removed_nodes$y_end,
                       # border = coloring_set$border_color,
                       border = rep("white", nrow(removed_nodes)),
                       # lty = coloring_set$lty_type,
                       lwd=2,
                       col = rep("white", nrow(removed_nodes))
        )
      }
      
    }
    
    #### continue from here
    if(edge_mapping) {
      
      if(custom_edge_mapping) {
        
        from <- graphnel_df$node_table[which(graphnel_df$node_table$node_id == highlight_nodes[1]),]
        to <- graphnel_df$node_table[which(graphnel_df$node_table$node_id == highlight_nodes[2]),]
        
        shape::Arrows(
          x0 = (from$x_start + from$x_end)/2,
          x1 = (to$x_start + to$x_end)/2,
          y0 = (from$y_start + from$y_end)/2,
          y1 = (to$y_start + to$y_end)/2,
          col = "red", 
          lwd=2, arr.length = 0.2, arr.type = "simple"
        )
        
      } else {
        if(!is.null(graphnel_df$edge_table)) {
          
          if(length(highlight_nodes) == 0) {
            shape::Arrows(
              x0 = graphnel_df$node_table[graphnel_df$edge_table$from,"x"],
              x1 = graphnel_df$node_table[graphnel_df$edge_table$to,"x"],
              y0 = graphnel_df$node_table[graphnel_df$edge_table$from,"y"],
              y1 = graphnel_df$node_table[graphnel_df$edge_table$to,"y"],
              col = graphnel_df$edge_table$color, 
              lty = ifelse(is.na(graphnel_df$edge_table$dashes), FALSE, graphnel_df$edge_table$dashes) + 1,
              lwd=2, arr.length = 0.2, arr.type = "simple"
            )
          } else {
            if(edge_in_mode) {
              index <- unique(which(graphnel_df$edge_table$from %in% highlight_nodes & graphnel_df$edge_table$to %in% highlight_nodes))
            } else {
              index <- unique(c(which(graphnel_df$edge_table$from %in% highlight_nodes), 
                                which(graphnel_df$edge_table$to %in% highlight_nodes)))
            }
            
            shape::Arrows(
              x0 = graphnel_df$node_table[graphnel_df$edge_table$from[index],"x"],
              x1 = graphnel_df$node_table[graphnel_df$edge_table$to[index],"x"],
              y0 = graphnel_df$node_table[graphnel_df$edge_table$from[index],"y"],
              y1 = graphnel_df$node_table[graphnel_df$edge_table$to[index],"y"],
              col = graphnel_df$edge_table$color[index], 
              lty = ifelse(is.na(graphnel_df$edge_table$dashes[index]), FALSE, graphnel_df$edge_table$dashes[index]) + 1,
              lwd=2, arr.length = 0.2, arr.type = "simple"
            )
          }
        }
      }
      
    }
    
    dev.off()
    return(img)
    
  }
  
  
  #### node color coding function ####
  node_color_generator <- function(pathway, sample_id = "mean", log_norm = TRUE, color_nodes = "psf_activities", custom_color_scale = NULL) {
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
        hover_text <- paste(paste("Log exp FC", round(log(pathway_exp_values + 0.00001), digits = 5)),
                            paste("Log PSF", round(log(pathway_psf_values + 0.00001), digits = 5)),
                            sep = "<br>")
        node_values <- cbind(round(log(pathway_exp_values + 0.00001), digits = 5), round(log(pathway_psf_values + 0.00001), digits = 5))
      } else {
        hover_text <- paste(paste("Exp FC", round(pathway_exp_values, digits = 5)),
                            paste("PSF", round(pathway_psf_values, digits = 5)),
                            sep = "<br>")
        node_values <- cbind(round(pathway_exp_values, digits = 5), round(pathway_psf_values, digits = 5))
      }
      
      names(hover_text) <- rownames(pathway$exp_fc)
      
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
    
    if(mapping_data) {
      if(color_nodes == "psf_activities") {
        pathway_node_values <- pathway_psf_values
      } else {
        pathway_node_values <- pathway_exp_values
      }
      
      if(is.null(custom_color_scale)) {
        node_colors <- psf:::color_code(values = pathway_node_values, pal1 = psf:::pal1, pal2 = psf:::pal2, log_scale = log_norm)
        
        node_colors <- data.frame(node_id = names(pathway_node_values)[c(which(pathway_node_values <= 0), which(pathway_node_values > 0))], 
                                  col = node_colors,
                                  text_col = unname(sapply(node_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                                  stringsAsFactors = F
        )
      } else {
        node_colors <- map_colors(values = pathway_node_values, colors = custom_color_scale)
        
        node_colors <- data.frame(node_id = names(pathway_node_values), 
                                  col = node_colors,
                                  text_col = unname(sapply(node_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                                  stringsAsFactors = F
        )
      }
      
      
      
      node_colors <- node_colors[which(node_colors$node_id %in% graph::nodes(pathway$graph)[which(unlist(graph::nodeData(pathway$graph, attr = "type")) != "map")]),]
      
    } else {
      stop("Please provide pathway with evalueated activity")
    }
    
    if(color_nodes == "psf_activities") {
      col_legend_title = ifelse(log_norm, "Log PSF value", "PSF value")
    } else {
      col_legend_title = ifelse(log_norm, "Log FC value", "FC value")
    }
    
    return(list(node_colors = node_colors, col_legend_title = col_legend_title, color_bar_lims = range(pathway_node_values), hover_text = hover_text, node_values = node_values))
  }
  
  map_colors <- function(values, colors) {
    middle_index <- length(colors) / 2
    negative_colors <- rev(colors[1:middle_index])
    positive_colors <- colors[(middle_index + 1):length(colors)]
    max_value <- max(abs(values))
    mapped_colors <- sapply(values, function(x) {
      if (x < 0) {
        # Map negative values to negative_colors
        index <- max(1, ceiling((abs(x) / max_value) * middle_index))
        return(negative_colors[index])
      } else if (x > 0) {
        # Map positive values to positive_colors
        index <- max(1, ceiling((x / max_value) * middle_index))
        return(positive_colors[index])
      } else {
        # Assign the middle color for zero values
        return(colors[middle_index])
      }
    })
    return(mapped_colors)
  }
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  #### visnet data extractor ####
  vis_extract <- function(vis_table, node_colors = NULL, higlight_node = NULL, custom_color_scale = NULL) {
    node_table <- vis_table$node_table[,c("node_id", "label", "vis_width", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title", "x", "y")]
    
    colnames(node_table) <- c("id", "label", "widthConstraint.maximum", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title", "x", "y")
    
    if(!is.null(node_colors)) {
      
      magick::autoviewer_disable() ### to avoid legend plotting before netowrk rendering
      legend_img <- magick::image_device(width = 480, height = 480)
      plot.new()
      
      if(is.null(custom_color_scale)) {
        psf:::color_legend_maker(x = 0.05, y = 0, leg = 0.9, cols = c(psf:::pal1(10), psf:::pal2(10)), title = node_colors$col_legend_title, lims = node_colors$color_bar_lims, digits=3, prompt=FALSE,
                                 lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
      } else {
        psf:::color_legend_maker(x = 0.05, y = 0, leg = 0.9, cols = c(custom_color_scale[1:50], custom_color_scale[51:100]), title = node_colors$col_legend_title, lims = node_colors$color_bar_lims, digits=3, prompt=FALSE,
                                 lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
      }
      
      temp_legend <- tempfile()
      
      magick::image_write(magick::image_trim(legend_img, fuzz = 0), path = temp_legend)
      
      legend_path <- paste('data:image/png;base64', RCurl::base64Encode(readBin(temp_legend, 'raw', file.info(temp_legend)[1, 'size']), 'txt'), sep = ',')
      
      legend_data_frame <- data.frame(
        id = as.character(max(as.integer(node_table$id)) + 1),
        label = "Color legend", 
        widthConstraint.maximum = NA, 
        font.color = "#000000", 
        size = 30, 
        font.size = 32, 
        borderWidth = 0, 
        color.border = "",
        color.background = "",
        shape = "image", 
        title = "",
        x = max(node_table$x),
        y = min(node_table$y) - 10,
        image = legend_path
      )
      
      node_table$image <- "unselected"
      node_table <- rbind(node_table, legend_data_frame)
      
      node_table$color.background <- unname(sapply(node_table$id, function(x) {
        if(x %in% node_colors$node_colors$node_id) {
          node_colors$node_colors[which(node_colors$node_colors$node_id == x),"col"]
        } else {
          "#BFFFBF"
        }
      }))
      
      node_table$color.border <- node_table$color.background
      if(!is.null(higlight_node)) {
        node_table$color.border[which(node_table$id %in% higlight_node)] <- "#48f542"
      }
      
      
      node_table$font.color <- unname(sapply(node_table$id, function(x) {
        if(x %in% node_colors$node_colors$node_id) {
          node_colors$node_colors[which(node_colors$node_colors$node_id == x),"text_col"]
        } else {
          "#000000"
        }
      }))
      
      node_table$title <- paste(node_table$title, node_colors$hover_text[node_table$id], sep = "<br>")
    }
    
    
    edge_table <- vis_table$edge_table[,c("id", "from", "to", "color", "arrows.to.enabled", "arrows.to.type", "label", "dashes")]
    
    if(all(is.na(node_table$x))) {
      node_table$x <- NULL
      node_table$y <- NULL
    }
    
    return(list(node_table = node_table, edge_table = edge_table))
  }
  
  if("image" %in% names(psf_saptial_results$spatial_psf_collection[[1]]$attrs)) {
    sink_name_to_id <- setNames(object = rownames(kegg_sink_to_process), nm = paste0(kegg_sink_to_process$Pathway_name, "->", kegg_sink_to_process$Sink)) 
  } else {
    sink_name_to_id <- unlist(lapply(names(psf_saptial_results$spatial_psf_collection), function(x) {
      setNames(object = paste0(psf_saptial_results$spatial_psf_collection[[x]]$sink.nodes, "; ", x), 
               nm = paste0(x, "->", unlist(graph::nodeData(psf_saptial_results$spatial_psf_collection[[x]]$graph, 
                                                           psf_saptial_results$spatial_psf_collection[[x]]$sink.nodes, attr = "label")))
      )
    }))
  }
  
  
  # Your existing code for data preparation
  # Extract the spatial image (stored as an array)
  spatial_image <- as.raster(psf_saptial_results$spatial_image)
  
  # Extract the low-resolution coordinates using GetTissueCoordinates
  coords <- psf_saptial_results$coords # GetTissueCoordinates(object = psf_spatial$spatial_psf_obj, scale = "lowres")
  
  coords$imagerow <- nrow(spatial_image) - coords$imagerow
  
  # Create a data frame for plotting
  spot_groups <- psf_saptial_results$meta.data$seurat_clusters # psf_spatial$spatial_psf_obj@meta.data$seurat_clusters
  # Define cluster colors
  unique_clusters <- sort(unique(spot_groups))
  cluster_colors <- setNames(gg_color_hue(length(unique_clusters)), unique_clusters)
  
  plotting_data <- data.frame(coords, spot_groups = spot_groups)
  
  # All points are highlighted initially
  plotting_data$Highlight <- "Highlighted"
  plotting_data$ColorCode <- cluster_colors[plotting_data$spot_groups]
  plotting_data$Opacity <- 1
  plotting_data$text_info <- paste(rownames(plotting_data), '\n', 'Cluster:', plotting_data$spot_groups)
  
  x_min_image <- 0
  x_max_image <- ncol(spatial_image)
  y_min_image <- 0
  y_max_image <- nrow(spatial_image)
  
  # Convert spot_groups to character for consistency
  spot_groups <- as.character(spot_groups)
  
  
  ui <- shinyUI(
    navbarPage("PSF Spatial Browser",
               tabPanel("App",
                        fluidPage(
                          tags$script('
        $(document).on("keydown", function (e) {
           if(e.which == 17) {
              Shiny.onInputChange("multiple_choice", "yes");
           }
        });
      '),
                          #### hover dialog design ####
                          tags$head(tags$style('
                           #hover_text {
                             position: absolute;
                             max-width: 600px;
                             min-width: 0px;
                             border-radius: 5px;
                             z-index: 100;
                             color: white;
                             background-color: rgba(0, 0, 0, 0.8);
                             }
                             ')),
                          #### hovering script ####
                          tags$script('
                  $(document).ready(function(){
                  // id of the plot
                  $("#pathway_image").mousemove(function(e){ 
                  
                  // ID of uiOutput
                  $("#hover_text").show();
                  $("#hover_text").css({
                  top: (e.pageY + 5) + "px",
                  left: (e.pageX + 5) + "px"
                  });     
                  });     
                  });
        '),
                          uiOutput("hover_text"),
                          # titlePanel(span("PSF spatial browser", style = "color:#1f992f"), windowTitle = "PSF SP"),
                          column(3,
                                 wellPanel(
                                   # selectInput(
                                   #   inputId = "annotation",
                                   #   label = "Select spot Annotation",
                                   #   choices = "",
                                   #   selected = "",
                                   #   selectize = FALSE,
                                   #   width = "100%"
                                   # ),
                                   selectInput(
                                     inputId = "cluster",
                                     label = "Select cluster",
                                     choices = "all", # c("all", as.character(unique_clusters)),
                                     selected = "all", # as.character(unique_clusters)[1],
                                     selectize = FALSE,
                                     width = "100%"
                                   ),
                                   tabsetPanel(
                                     id = "app_mode",
                                     tabPanel(
                                       title = "Top features",
                                       selectizeInput(
                                         inputId = "selected_feature",
                                         label = "Feature to visualize",
                                         choices = NULL,
                                         selected = NULL,
                                         width = "100%"
                                       ),
                                       htmlOutput('feature_stat')
                                     ),
                                     tabPanel(
                                       title = "Pathway vis",
                                       selectizeInput("selected_pathway", label = "Select pathway", choices = NULL, options = list(
                                         placeholder = 'Please select a pathway',
                                         onInitialize = I('function() { this.setValue(""); }')
                                       )),
                                       radioButtons(
                                         "node_coloring_type",
                                         label = "Select mapping value",
                                         choices = list("PSF" = 1, "FC" = 2),
                                         selected = 1,
                                         inline = TRUE
                                       )
                                     )
                                   ),
                                   htmlOutput('error_text'),
                                   width = "100%"
                                 )
                          ),
                          column(
                            9,
                            plotlyOutput("spatial_plot", height = "550px", width = "650px"),
                            tabsetPanel(id = "net_vis_panels",
                                        tabPanel(title='KEGG', value = 1, useShinyjs(),
                                                 imageOutput(
                                                   "pathway_image",
                                                   dblclick = "image_click",
                                                   brush = brushOpts(id = "image_brush", delayType = "debounce", resetOnNew = TRUE),
                                                   hover = hoverOpts(id = "image_hover", delay = 0)
                                                 )        
                                        ),
                                        tabPanel(title='VisNet',
                                                 visNetworkOutput("visnet", height = "800px")
                                        )
                            )
                          )
                        )
               ),
               tabPanel("Help",
                        tabsetPanel(
                          tabPanel("Welcome", includeMarkdown(system.file('help_pages/welcome.md', package='psf'))),
                          # tabPanel("Getting Started", includeMarkdown("help_pages/getting_started.md")),
                          tabPanel("User Interface Guide", includeMarkdown(system.file('help_pages/ui_guide.md', package='psf'))),
                          # tabPanel("Features", includeMarkdown("help_pages/features.md")),
                          tabPanel("Use Case", includeMarkdown(system.file('help_pages/use_case.md', package='psf'))),
                          tabPanel("FAQ & Troubleshooting", includeMarkdown(system.file('help_pages/faq.md', package='psf'))),
                          tabPanel("Resources", includeMarkdown(system.file('help_pages/resources.md', package='psf'))),
                          tabPanel("Contact", includeMarkdown(system.file('help_pages/contact.md', package='psf') ))
                        )
               )
    )
  )
  
  
  server <- function(input, output, session) {
    # Reactive values
    v <- reactiveValues(
      feature = NULL,
      plotly_data = plotting_data,
      plotly_title = "PSF clusters",
      value_type = "psf",
      spatial_plotly = NULL,
      graphnel_df = NULL,
      pathway_image = NULL,
      pathway_name = NULL,
      node_selection_error = "",
      plotly_selection = NULL,
      selected_feature = "",
      highlight_feature_id = NULL,
      vis_samples = NULL,
      feature_metdata = NULL,
      selected_node_name = setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
    )
    
    updateSelectInput(session, "cluster", choices = c("all", as.character(unique_clusters)), selected = "all")
    
    if(!("image" %in% names(psf_saptial_results$spatial_psf_collection[[1]]$attrs))) {
      hideTab(session, inputId = "net_vis_panels", target = "1")
    }
    
    updateSelectizeInput(session, "selected_pathway", choices = names(psf_saptial_results$spatial_psf_collection), server = FALSE)
    
    #### rendering spatial plotly ####
    observe({
      # Create the Plotly plot
      v$spatial_plotly <- plot_ly(
        data = v$plotly_data,
        x = ~imagecol,
        y = ~imagerow,
        customdata = rownames(v$plotly_data),
        type = 'scatter',
        mode = 'markers',
        text = ~text_info,
        hoverinfo = 'text',
        marker = list(
          size = 6,
          color = ~ColorCode,
          opacity = ~Opacity
        )
      ) %>%
        layout(
          title = list(text = v$plotly_title, x = 0.2),
          images = list(
            list(
              source = raster2uri(spatial_image),
              xref = "x",
              yref = "y",
              x = x_min_image,
              y = y_max_image,
              sizex = (x_max_image - x_min_image),
              sizey = (y_max_image - y_min_image),
              sizing = "stretch",
              opacity = 1,
              layer = "below"
            )
          ),
          xaxis = list(
            title = "",
            range = c(x_min_image, x_max_image),
            scaleanchor = "y",
            showgrid = FALSE,
            zeroline = FALSE,
            showticklabels = FALSE
          ),
          yaxis = list(
            title = "",
            range = c(y_min_image, y_max_image),
            scaleanchor = "x",
            showgrid = FALSE,
            zeroline = FALSE,
            showticklabels = FALSE
          ),
          showlegend = FALSE,  # Hide legend to prevent clutter
          dragmode = "pan"
        ) %>%
        config(scrollZoom = TRUE)
    })
    
    
    # Observe the selected cluster and update the plot accordingly
    observeEvent(input$cluster, {
      
      v$selected_feature <- ""
      v$feature_metdata <- NULL
      v$plotly_selection <- NULL
      v$highlight_feature_id <- NULL
      v$plotly_title <- paste0("PSF cluster ", input$cluster)
      
      if(input$cluster != "all") {
        updateSelectizeInput(inputId = "selected_feature", 
                             choices = c("", paste0(psf_saptial_results$psf_ident_markers[[input$cluster]]$Pathway_name, "->", psf_saptial_results$psf_ident_markers[[input$cluster]]$Sink)), 
                             selected = "", server = TRUE)
      } else {
        updateSelectizeInput(inputId = "selected_feature", 
                             choices = NULL, server = TRUE)
      }
      
      if(input$cluster == "all") {
        plotting_data$ColorCode <- cluster_colors[plotting_data$spot_groups]
        plotting_data$Opacity <- 1
      } else {
        # Add a column to indicate whether to highlight the point
        plotting_data$Highlight <- ifelse(
          plotting_data$spot_groups == input$cluster, "Highlighted", "Normal"
        )
        
        # Set color and opacity based on whether the point is highlighted
        plotting_data$ColorCode <- ifelse(
          plotting_data$Highlight == "Highlighted",
          cluster_colors[plotting_data$spot_groups],
          "grey"
        )
        
        plotting_data$Opacity <- ifelse(
          plotting_data$Highlight == "Highlighted",
          1,
          0.3
        )
      }
      
      plotting_data$text_info <- paste(rownames(plotting_data), '\n', 'Cluster:', plotting_data$spot_groups)
      
      v$plotly_data <- plotting_data
      
    })
    
    observeEvent(input$selected_feature, {
      v$selected_feature <- input$selected_feature
    })
    
    # Initial spatial plot rendering
    output$spatial_plot <- renderPlotly({
      v$spatial_plotly %>%
        event_register()
    })
    
    #### plotly selection events ####
    observeEvent(event_data("plotly_selected"), {
      v$plotly_selection <- event_data("plotly_selected")$customdata
    })
    
    observeEvent(event_data("plotly_deselect"), {
      v$plotly_selection <- event_data("plotly_selected")$customdata
    })
    
    #### pathway plotting variables observer ####
    observeEvent(c(v$selected_feature, input$selected_pathway, input$app_mode, v$plotly_selection, input$cluster), {
      
      if(!is.null(v$plotly_selection)) {
        v$vis_samples <- v$plotly_selection
      } else {
        if(input$cluster == "all") {
          v$vis_samples <- rownames(psf_saptial_results$meta.data)
        } else {
          v$vis_samples <- rownames(psf_saptial_results$meta.data)[which(as.character(psf_saptial_results$meta.data$seurat_clusters) == input$cluster)]
        }
      }
      
      if(input$app_mode == "Top features") {
        if(v$selected_feature != "") {
          v$highlight_feature_id <- unlist(strsplit(sink_name_to_id[v$selected_feature], split = "; "))[1]
          v$pathway_name <- unlist(strsplit(sink_name_to_id[v$selected_feature], split = "; "))[2]
          v$plotly_title <- paste0(v$selected_feature, " PSF")
        }
        
      } else {
        if(!is.null(input$selected_pathway)) {
          v$pathway_name <- input$selected_pathway
        }
      }
      
    })
    
    observeEvent(v$selected_feature, {
      if(input$cluster != "all" & v$selected_feature != "") {
        metdata_text <- psf_saptial_results$psf_ident_markers[[input$cluster]][gsub("_", "-", sink_name_to_id[v$selected_feature]),c("avg_log2FC", "p_val_adj", "Process", "Dowstream_pathway")]
        
        v$feature_metdata <- paste0(
          paste("<font color=\"#000000\"><b>",
                "Avg log2FC: ", round(metdata_text$avg_log2FC, digits = 3),
                "; P_adj: ", metdata_text$p_val_adj, "<br>",
                "</b></font>"), 
          paste("<font color=\"#1f992f\"><b>",
                "Process: ", metdata_text$Process,  "<br>",
                "Downstream pathway:", metdata_text$Dowstream_pathway,
                "</b></font>")
        )
      } else {
        v$feature_metdata <- NULL
      }
    })
    
    #### Selected cluster feature information #### 
    output$feature_stat <- renderText({
      v$feature_metdata
    })
    
    
    #### colored featureplot data generator ####
    observeEvent(v$highlight_feature_id, {
      feature_values <- log(psf_saptial_results$spatial_psf_collection[[v$pathway_name]][["psf_activities"]][v$highlight_feature_id,])
      
      plotting_data <- data.frame(coords, spot_groups = spot_groups)
      plotting_data$ColorCode <- map_colors(values = feature_values, colors = Seurat:::SpatialColors(n = 100)) # Seurat:::SpatialColors(n = 100)[cut(feature_values, breaks = 100, labels = FALSE)]
      plotting_data$Opacity <- 1
      plotting_data$text_info <- paste(rownames(plotting_data), '\n', 'Cluster:', plotting_data$spot_groups, '\n', "Log PSF", round(feature_values, digits = 5))
      v$plotly_data <- plotting_data
    })
    
    #### pathway coloring and plotting observer ####
    observeEvent(c(v$pathway_name, v$vis_samples, v$highlight_feature_id, v$value_type, input$node_coloring_type), {
      if(!is.null(v$pathway_name)) {
        if(v$pathway_name != "") {
          clust_pathway <- psf_subset(psf_collection = psf_saptial_results$spatial_psf_collection, 
                                      pathway_name = v$pathway_name, 
                                      sample_list = v$vis_samples)
          
          if(!is.null(clust_pathway)) {
            v$node_colors <- node_color_generator(pathway = clust_pathway, log_norm = T, sample_id = "mean", 
                                                  color_nodes = c("psf_activities", "fold_change")[as.integer(input$node_coloring_type)], 
                                                  custom_color_scale = Seurat:::SpatialColors(n = 100))
          }
          
          v$graphnel_df <- graphnel_to_df(clust_pathway, extended = TRUE)
          
          v$graphnel_df$node_table[rownames(v$node_colors$node_values), "expression"] <- v$node_colors$node_values[,1]
          v$graphnel_df$node_table[rownames(v$node_colors$node_values), "signal"] <- v$node_colors$node_values[,2]
          
          
          if("image" %in% names(psf_saptial_results$spatial_psf_collection[[1]]$attrs)) {
            img_path <- system.file("extdata", "old_imgs", paste0(gsub("path:", "", clust_pathway$attrs$name), ".png"), package = "psf")
            
            v$pathway_image <- magick::image_read(img_path)
            
            v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = NULL, pathway_image = v$pathway_image,
                                              node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, 
                                              highlight_nodes = v$highlight_feature_id, highlight_color = "green", 
                                              col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, 
                                              present_node_modifications = FALSE, removed_nodes = NULL, custom_color_scale = Seurat:::SpatialColors(n = 100)) %>% 
              image_write(tempfile(fileext='png'), format = 'png')
          }
          
        }
      }
    })
    
    
    #### image click events ####
    observeEvent(input$image_click, {
      
      if(!is.null(input$multiple_choice)) {
        if(input$multiple_choice == "yes") {
          v$selected_node_name <- unique(rbind(v$selected_node_name, v$graphnel_df$node_table[which(input$image_click[[1]] <= v$graphnel_df$node_table$x_end & input$image_click[[1]] >= v$graphnel_df$node_table$x_start &
                                                                                                      input$image_click[[2]] <= v$graphnel_df$node_table$y_end & input$image_click[[2]] >= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]))
        } else {
          v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
          v$selected_node_name <- v$graphnel_df$node_table[which(input$image_click[[1]] <= v$graphnel_df$node_table$x_end & input$image_click[[1]] >= v$graphnel_df$node_table$x_start &
                                                                   input$image_click[[2]] <= v$graphnel_df$node_table$y_end & input$image_click[[2]] >= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]
        }
      } else {
        v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
        v$selected_node_name <- v$graphnel_df$node_table[which(input$image_click[[1]] <= v$graphnel_df$node_table$x_end & input$image_click[[1]] >= v$graphnel_df$node_table$x_start &
                                                                 input$image_click[[2]] <= v$graphnel_df$node_table$y_end & input$image_click[[2]] >= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]
      }
      
      
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = NULL, pathway_image = v$pathway_image,
                                        node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                        col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = FALSE, removed_nodes = NULL, custom_color_scale = Seurat:::SpatialColors(n = 100)) %>%
        image_write(tempfile(fileext='png'), format = 'png')
      
    })
    
    #### image brush events ####
    observeEvent(input$image_brush, {
      if(!is.null(input$multiple_choice)) {
        if(input$multiple_choice == "yes") {
          
          v$selected_node_name <- unique(rbind(v$selected_node_name, v$graphnel_df$node_table[which(input$image_brush[1:4]$xmax >= v$graphnel_df$node_table$x_end & input$image_brush[1:4]$xmin <= v$graphnel_df$node_table$x_start &
                                                                                                      input$image_brush[1:4]$ymax >= v$graphnel_df$node_table$y_end & input$image_brush[1:4]$ymin <= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]))
          
        } else {
          v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
          v$selected_node_name <- v$graphnel_df$node_table[which(input$image_brush[1:4]$xmax >= v$graphnel_df$node_table$x_end & input$image_brush[1:4]$xmin <= v$graphnel_df$node_table$x_start &
                                                                   input$image_brush[1:4]$ymax >= v$graphnel_df$node_table$y_end & input$image_brush[1:4]$ymin <= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]
        }
      } else {
        v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
        v$selected_node_name <- v$graphnel_df$node_table[which(input$image_brush[1:4]$xmax >= v$graphnel_df$node_table$x_end & input$image_brush[1:4]$xmin <= v$graphnel_df$node_table$x_start &
                                                                 input$image_brush[1:4]$ymax >= v$graphnel_df$node_table$y_end & input$image_brush[1:4]$ymin <= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]
      }
      
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = NULL, pathway_image = v$pathway_image,
                                        node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                        col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = FALSE, removed_nodes = NULL, custom_color_scale = Seurat:::SpatialColors(n = 100)) %>%
        image_write(tempfile(fileext='png'), format = 'png')
      
      
    })
    
    observeEvent(input$clicked_node$nodes, {
      v$selected_node_name <- v$graphnel_df$node_table[which(v$graphnel_df$node_table$node_id == unlist(input$clicked_node$nodes)),c("node_id", "component_id_s", "label")]
    })
    
    
    #### plotly visualization by node selection ####
    observeEvent(c(v$selected_node_name, input$node_coloring_type), {
      if(nrow(v$selected_node_name) == 1) {
        
        plotting_data <- data.frame(coords, spot_groups = spot_groups)
        
        feature_values <- log(psf_saptial_results$spatial_psf_collection[[v$pathway_name]][[c("psf_activities", "exp_fc")[as.integer(input$node_coloring_type)]]][v$selected_node_name$node_id,])
        
        plotting_data$ColorCode <-  map_colors(values = feature_values, colors = Seurat:::SpatialColors(n = 100)) # Seurat:::SpatialColors(n = 100)[cut(feature_values, breaks = 100, labels = FALSE)]
        plotting_data$Opacity <- 1
        plotting_data$text_info <- paste(rownames(plotting_data), '\n', 'Cluster:', plotting_data$spot_groups, '\n', "Log PSF", round(feature_values, digits = 5))
        
        v$plotly_data <- plotting_data
        
        v$plotly_title <- paste0(v$selected_node_name$label, c(" Log PSF", " Log Exp FC")[as.integer(input$node_coloring_type)])
      } else {
        if(nrow(v$selected_node_name) > 1) {
          v$node_selection_error <- "Only one node can be visualized"
          delay(3000, v$node_selection_error <- "")
        }
      }
    })
    
    output$error_text <- renderText({
      paste("<font color=\"#ff0000\"><b>", v$node_selection_error, "</b></font>")
    })
    
    #### image hover events ####
    output$hover_text <- renderUI({
      h5(textOutput("hover_data"))
    })
    
    output$hover_data <- renderText({
      if(length(which(input$image_hover[[1]] <= v$graphnel_df$node_table$x_end & input$image_hover[[1]] >= v$graphnel_df$node_table$x_start &
                      input$image_hover[[2]] <= v$graphnel_df$node_table$y_end & input$image_hover[[2]] >= v$graphnel_df$node_table$y_start)) > 0 ) {
        
        gr_name <- v$graphnel_df$node_table[which(input$image_hover[[1]] <= v$graphnel_df$node_table$x_end & input$image_hover[[1]] >= v$graphnel_df$node_table$x_start &
                                                    input$image_hover[[2]] <= v$graphnel_df$node_table$y_end & input$image_hover[[2]] >= v$graphnel_df$node_table$y_start),"title"]
        
        if(!is.null(v$node_colors)) {
          node_id <- v$graphnel_df$node_table[which(input$image_hover[[1]] <= v$graphnel_df$node_table$x_end & input$image_hover[[1]] >= v$graphnel_df$node_table$x_start &
                                                      input$image_hover[[2]] <= v$graphnel_df$node_table$y_end & input$image_hover[[2]] >= v$graphnel_df$node_table$y_start),"node_id"]
          gr_name <- paste0(gr_name, "\n", v$node_colors$hover_text[node_id])
        }
        
        ### fix next line issue
        gsub("<br>", "\n", gr_name)
        
      }
    })
    
    
    output$pathway_image <- renderImage({
      return(list(src = ifelse(is.null(v$graphnel_df), "", v$image_file), contentType = "image/png", alt = "pathway_image"))
    }, deleteFile = TRUE)
    
    
    #### visnet rendering ####
    output$visnet <- renderVisNetwork({
      if(!is.null(v$graphnel_df)) {
        visNetwork(nodes = vis_extract(v$graphnel_df, node_colors = v$node_colors, higlight_node = v$highlight_feature_id, custom_color_scale = Seurat:::SpatialColors(n = 100))$node_table, 
                   edges = vis_extract(v$graphnel_df)$edge_table, width = "100%", height = "800px") %>% 
          visIgraphLayout(layout = "layout_nicely") %>% 
          visInteraction(navigationButtons = TRUE, multiselect = F, keyboard = F, selectConnectedEdges = T) %>%
          visOptions(highlightNearest = T) %>%
          visExport() %>%
          visEvents(click = "function(nodes) {
                      console.info('click')
                      console.info(nodes)
                      Shiny.onInputChange('clicked_node', {nodes : nodes.nodes, edges : nodes.edges});
                      ;}"
          )
      }
      
    })
    
  }
  
  
  shinyApp(ui = ui, server = server)
  
}

