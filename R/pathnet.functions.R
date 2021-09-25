#' Map gene data onto a pathway graph
#'
#'
#' @param g Tha pathway graph of graphNEL class
#' @param entrez.fc gene expression fold change matrix with entrez gene rownames (A matrix data associated with the genes; rownames represent genes; a single gene-row may contain one or many data values;)
#'
#' @return graphNEL object with the column-wise averaged gene data kept in the nodedata attribute of the graph
map.gene.data <- function(g, entrez.fc){
  gene.data <- graph::nodeData(g)
  gene.data <- lapply(gene.data, function (x,y) {
    genes.in.node = which(rownames(y) %in% x$genes)
    expression.values <- y[genes.in.node,]
    #         show(length(expression.values))
    if (length(expression.values)>0){
      if (x$type == "gene"){
        expression <- mean(expression.values)
        #         show(expression.values)
      }
      else { if(x$type == "group")
        expression <- min(expression.values)
      #         show(expression.values)
      }
    }
    else expression <-1

    #     show(expression)
    x$expression <- expression
    return(x)
  }, entrez.fc)

  for(i in 1:length(gene.data)){
    if(gene.data[[i]]$type == "gene")
      graph::nodeData(g, names(gene.data)[i], "expression") =  gene.data[[i]]$expression
  }
  return(g)
}

ensembl.id.conversion <- function(indata.FC, convert.from = "ensembl_gene_id", convert.to = "entrezgene"){

  z = rownames(indata.FC)

  cat("Loading biomart data table\n")
  #mart<-biomaRt::useMart('ensembl')
  mart = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  mart<-biomaRt::useDataset("hsapiens_gene_ensembl", mart=mart)
  biomart.table = biomaRt::getBM( c(convert.from, convert.to), values = z, mart = mart,
                                  filters=convert.from, uniqueRows = T )

  na.ind = which(apply(apply(biomart.table, 2, is.na),1,any))
  if(length(na.ind) > 0)
    biomart.table = biomart.table[-na.ind, ]

  Entrez.FC = matrix(0, nrow = nrow(biomart.table), ncol = ncol(indata.FC))
  colnames(Entrez.FC) <- colnames(indata.FC)
  rownames(Entrez.FC) <- biomart.table[,2]

  cat("Performing ID conversion:\n")
  percent5 = base::floor(nrow(biomart.table)/20)
  # cat(percent5)
  perc = 0
  for(i in 1:nrow(biomart.table)) {
    if((i %% percent5) == 0){
      perc = perc + 5
      # show(i)
      cat(perc, "%\n")
    }
    gene.ID = biomart.table[i,1]
    #     Entrez.FC[i,] = indata.FC[which(tolower(rownames(indata.FC)) == tolower(gene.ID)),]
    row.index = which(tolower(z) == tolower(gene.ID))
    # if(length(row.index) != 0){
    if(length(row.index) > 1 )
      Entrez.FC[i,] = colMeans(indata.FC[row.index,])
    else
      Entrez.FC[i,] = indata.FC[row.index,]
    # }
  }


  return(Entrez.FC)
}

psf.from.env.indata <- function(indata.fc){
  entrez.fc = ensembl.id.conversion(indata.fc)
  return(psf.env.entrez.fc(entrez.fc))
}

psf.from.env.entrez.fc <- function(entrez.fc, kegg.collection, split = TRUE, calculate.significance = T, bst.steps = 200, sum = FALSE){
  psf.results.collection = list()
  psf.results.collections = list()
  for(c in 1:ncol(entrez.fc)){
    # cat("Column", c)
    for(i in 1:length(kegg.collection)){
      cat("\nKEGG collection number: ", i, " exp matrix column: ", c)
      pathway = names(kegg.collection)[i]
      # show(pathway)
      if(!length(kegg.collection[[i]])==0){
        # cat("PSF computation on: ", pathway, "\n")
        entrez.column = as.matrix(entrez.fc[,c])
        g = map.gene.data(kegg.collection[[i]]$graph, entrez.column)

        psf.results.collection[[pathway]] = psf.flow(g, kegg.collection[[i]]$order, kegg.collection[[i]]$sink.nodes, split, sum = sum)

        psf.results.collection[[pathway]] = c(psf.results.collection[[pathway]],
                                              attrs=list(kegg.collection[[i]]$attrs),
                                              order=list(kegg.collection[[i]]$order))

      }
    }
    if(calculate.significance){
      cat("\nPerforming bootstrap calculations with", bst.steps, " steps\n")
      psf.results = bootstrap.significance(psf.results.collection, entrez.fc, bst.steps)
      psf.results.processed = process.psf.results(psf.results)
    } else {
      psf.results.processed = psf.results.collection
    }


    psf.results.collections[[c]] = psf.results.processed
  }

  return(psf.results.collections)
}

bootstrap.significance <- function(psf.results, entrez.fc, bst.steps){
  for(i in 1:length(psf.results)){
    result.it = psf.results[[i]]
    cat("Performing bootstrap on: ", result.it$attrs$title, "\n")
    boot.mat = psf.flow.boot(result.it$graph, entrez.fc,
                             result.it$order, result.it$eval.exprs,
                             result.it$sink.nodes, result.it$I, bst.steps)
    psf.results[[i]] = c(psf.results[[i]], list(boot.mat=boot.mat))
  }
  return(psf.results)
}


psf.from.env.to.table <- function(env){
  if(is.null(env))
    stop("null env argument")

  if(!(kegg.collection %in% ls())){
    cat("loading kegg collection\n")
    load("data_to_be_deleted/kegg.collection/kegg.collection.Rdata")
  }
  if(length(kegg.collection) == 0)
    stop("Empty kegg.collection")
  cat("kegg.collection with ", length(kegg.collection), " pathways successfully loaded\n")

  cat("PSF computation started.\n")

  #genewise normalized data: considered as logFC.
  indata = env$indata

  #anti-log the indata to get the FC values
  indata.fc = 10^indata

  psf.results.collection = psf.from.env.indata(indata.fc)

}

process.psf.results <- function(psf.results){
  for(i in 1:length(psf.results)){
    if(length(psf.results[[i]]$signal.at.sink)>0){
      mean.psf = base::mean(unlist(psf.results[[i]]$signal.at.sink))
      psf.results[[i]]$mean.psf = mean.psf
      max.psf = base::max(unlist(psf.results[[i]]$signal.at.sink))
      psf.results[[i]]$max.psf = max.psf
      max.sinks = which(psf.results[[i]]$signal.at.sink == max.psf)
      psf.results[[i]]$max.sinks = names(max.sinks)
      if(!is.null(psf.results[[i]]$boot.mat)){
        p.values = sig.calc(psf.results[[i]]$signal.at.sink, psf.results[[i]]$boot.mat)
        psf.results[[i]]$p.values = p.values
        psf.results[[i]]$max.p.values = p.values[max.sinks]
        psf.results[[i]]$mean.p.value = 10^(base::mean(log10(p.values)))
      }

    } else {
      # cat("no signals at ", psf.results[[i]]$attrs$title)
    }
  }
  return(psf.results)
}


determine.sink.nodes.for.collection <- function(kegg.collection){
  for(i in 1:length(kegg.collection)){
    pathway = kegg.collection[[i]]
    cat(pathway$attrs$title, "\n")
    sink.nodes = determine.sink.nodes(pathway)
    kegg.collection[[i]]$sink.nodes = sink.nodes
  }
  return(kegg.collection)
}

determine.sink.nodes <- function(pathway){
  sink.nodes =  NULL
  g = pathway$graph
  order = pathway$order$node.order
  for(node in graph::nodes(g)){
    node.children = graph::edges(g)[[node]]
    parent.nodes = NULL
    for(parent in graph::nodes(g)){
      if(node %in% graph::edges(g)[[parent]])
        parent.nodes = c(parent, parent.nodes)
    }
    isSinkNode = F
    if(length(parent.nodes) != 0){
      if(length(node.children) == 0)
        isSinkNode = T
      #         else{
      #           node.rank = order[[node]]
      #           for(child in node.children){
      #             isSinkNode = T
      #             child.rank = order[[child]]
      #             if(child.rank <= node.rank)
      #               isSinkNode = F
      #           }
      #         }
    }
    if(isSinkNode)
      sink.nodes = c(sink.nodes, node)
  }

  return(sink.nodes)
}



