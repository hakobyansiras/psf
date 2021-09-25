plot.pathway.henry <- function(graph, sink.nodes = NULL){
  g <- igraph::igraph.from.graphNEL(graph)

  nodedata <- graph@nodeData@data
  edgedata <- graph@edgeData@data
  node.type <- sapply( nodedata, function(x) x$kegg.type )
  edge.subtype <- sapply( edgedata, function(x) c( x$subtype1,x$subtype2) )

  edge.subtype <- lapply( edge.subtype, function(x) if("phosphorylation"%in%x && sum(!is.na(x))==1) c(x,"actiigraph::Vation") else x )
  edge.subtype <- lapply( edge.subtype, function(x) if("dephosphorylation"%in%x && sum(!is.na(x))==1) c(x,"actiigraph::Vation") else x )
  edge.subtype <- lapply( edge.subtype, function(x) if("ubiquitination"%in%x && sum(!is.na(x))==1) c(x,"actiigraph::Vation") else x )
  edge.subtype <- lapply( edge.subtype, function(x) if("methylation"%in%x && sum(!is.na(x))==1) c(x,"actiigraph::Vation") else x )

  title.node <- grep( "TITLigraph::E:", igraph::V(g)$label )
  if("sink.nodes" %in% ls())
    sink.nodes <- match( sink.nodes, names(nodedata) )
  gene.nodes <- which(node.type%in%c("gene","ortholog"))
  comp.nodes <- which(node.type%in%c("compound"))
  spare.nodes <- setdiff( seq(igraph::V(g)), c(title.node,gene.nodes,comp.nodes) )


  igraph::V(g)$names <- igraph::V(g)$label

  igraph::V(g)$x <- sapply( nodedata, function(x) as.numeric( x$kegg.gr.x ) )
  igraph::V(g)$y <- -sapply( nodedata, function(x) as.numeric( x$kegg.gr.y ) )

  igraph::V(g)$shape <- "rectangle"
  igraph::V(g)$shape[comp.nodes] <- "circle"

  igraph::V(g)$frame.color <- "black"
  igraph::V(g)$frame.color[spare.nodes] <- "gray60"

  igraph::V(g)$label.color <- "black"
  igraph::V(g)$label.color[spare.nodes] <- "gray60"

  igraph::V(g)$label[title.node] <- sub( "TITLigraph::E:", "", igraph::V(g)$label[title.node] )
  igraph::V(g)$label <- sapply( igraph::V(g)$label, function(x)
  {
    if( nchar(x) > 25 )
    {
      s <- strsplit(x," ")[[1]]
      return( paste( paste( s[ 1:ceiling(length(s)/2) ], collapse=" "  ),
                     paste( s[ (ceiling(length(s)/2)+1):length(s) ], collapse=" "  ), sep="\n" ) )

    }else return( x )
  })
  igraph::V(g)$label[title.node] = paste( igraph::V(g)$label[title.node], sep="\n")

  igraph::V(g)$label.cex <- 0.6
  igraph::V(g)$label.cex[title.node] <- 1.2
  igraph::V(g)$label.cex[spare.nodes] <- 0.6

  igraph::V(g)$color <- "white"
  igraph::V(g)$color[title.node] <- "gray75"
  #igraph::V(g)$color[gene.nodes] <- colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000)[999*(signal.Values[gene.nodes]-signal.Values.lim[1])/(signal.Values.lim[2]-signal.igraph::Values.lim[1])+1]

  igraph::V(g)$size <- 48
  igraph::V(g)$size2 <- 20
  igraph::V(g)$size[gene.nodes] <- 10
  igraph::V(g)$size2[gene.nodes] <- 6
  igraph::V(g)$size[spare.nodes] <- 28
  igraph::V(g)$size2[spare.nodes] <- 8
  igraph::V(g)$size[comp.nodes] <- 8

  igraph::E(g)$color <- "gray20"
  igraph::E(g)$label.color <- "black"
  igraph::E(g)$arrow.size <- 0
  igraph::E(g)$arrow.width <- 0
  igraph::E(g)$lty <- 3

  igraph::E(g)$label <- ""
  igraph::E(g)$label[which( sapply( edge.subtype[apply( igraph::get.edgelist(g), 1, paste, collapse="|" )], function(x) "phosphorylation" %in% x ) )] <- "+p"
  igraph::E(g)$label[which( sapply( edge.subtype[apply( igraph::get.edgelist(g), 1, paste, collapse="|" )], function(x) "dephosphorylation" %in% x ) )] <- "-p"
  igraph::E(g)$label[which( sapply( edge.subtype[apply( igraph::get.edgelist(g), 1, paste, collapse="|" )], function(x) "ubiquitination" %in% x ) )] <- "+u"
  igraph::E(g)$label[which( sapply( edge.subtype[apply( igraph::get.edgelist(g), 1, paste, collapse="|" )], function(x) "methylation" %in% x ) )] <- "+m"
  igraph::E(g)$label[which( sapply( edge.subtype[apply( igraph::get.edgelist(g), 1, paste, collapse="|" )], function(x) "dissociation" %in% x ) )] <- "#"


  if("highlight.genes" %in% ls() && !is.null(highlight.genes))
  {
    igraph::V(g)$label[which(highlight.genes!="")] <- paste(igraph::V(g)[which(highlight.genes!="")]$label," (",highlight.genes[which(highlight.genes!="")],")",sep="")
    igraph::V(g)$frame.color[which(highlight.genes!="")] <- "indianred1"
  }
  if("highlight.sinks" %in% ls())
    if(highlight.sinks) igraph::V(g)$color[sink.nodes] <- "lightblue2"

  plot(g,asp=0)


  igraph::V(g)$label[spare.nodes] <- ""
  igraph::V(g)$color[spare.nodes] <- NA
  igraph::V(g)$frame.color[spare.nodes] <- NA


  # --> edges
  g2 <- g - igraph::E(g)[ which( !sapply( edge.subtype[apply( igraph::get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("actiigraph::Vation","expression","reaction") %in% x) ) )  ]
  igraph::E(g2)$arrow.size <- 0.4
  igraph::E(g2)$arrow.width <- 0.8
  igraph::E(g2)$lty <- 1

  par(new=T)
  plot(g2,asp=0)



  # ..> edges
  g2 <- g - igraph::E(g)[ which( !sapply( edge.subtype[apply( igraph::get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("indirect effect") %in% x) ) )  ]
  igraph::E(g2)$arrow.size <- 0.4
  igraph::E(g2)$arrow.width <- 0.8
  igraph::E(g2)$lty <- 2

  par(new=T)
  plot(g2,asp=0)



  # --| edges
  g2 <- g - igraph::E(g)[ which( !sapply( edge.subtype[apply( igraph::get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("inhibition","repression") %in% x) ) )  ]
  igraph::E(g2)$arrow.size <- 0.1
  igraph::E(g2)$arrow.width <- 16
  igraph::E(g2)$lty <- 1

  par(new=T)
  plot(g2,asp=0)



  # -- edges
  g2 <- g - igraph::E(g)[ which( !sapply( edge.subtype[apply( igraph::get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("binding/association","dissociation","compound") %in% x) ) )  ]
  igraph::E(g2)$arrow.size <- 0
  igraph::E(g2)$arrow.width <- 0
  igraph::E(g2)$lty <- 1

  par(new=T)
  plot(g2,asp=0)
  igraph::tkplot(g2)


}
