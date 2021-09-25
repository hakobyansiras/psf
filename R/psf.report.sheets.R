psf.signal.sheets <- function(psf.results)
{

  plot.psf.pathway <- function( psf.object, signal.values, signal.values.lim, main="",
                                highlight.sinks=FALSE, highlight.genes=NULL )
  {
    g <- igraph.from.graphNEL(psf.object$graph)

    nodedata <- psf.object$graph@nodeData@data
    edgedata <- psf.object$graph@edgeData@data
    node.type <- sapply( nodedata, function(x) x$kegg.type )
    edge.subtype <- sapply( edgedata, function(x) c( x$subtype1,x$subtype2) )

    edge.subtype <- lapply( edge.subtype, function(x) if("phosphorylation"%in%x && sum(!is.na(x))==1) c(x,"activation") else x )
    edge.subtype <- lapply( edge.subtype, function(x) if("dephosphorylation"%in%x && sum(!is.na(x))==1) c(x,"activation") else x )
    edge.subtype <- lapply( edge.subtype, function(x) if("ubiquitination"%in%x && sum(!is.na(x))==1) c(x,"activation") else x )
    edge.subtype <- lapply( edge.subtype, function(x) if("methylation"%in%x && sum(!is.na(x))==1) c(x,"activation") else x )

    title.node <- grep( "TITLE:", V(g)$label )
    sink.nodes <- match( psf.object$sink.nodes, names(nodedata) )
    gene.nodes <- which(node.type%in%c("gene","ortholog"))
    comp.nodes <- which(node.type%in%c("compound"))
    spare.nodes <- setdiff( seq(V(g)), c(title.node,gene.nodes,comp.nodes) )


    V(g)$names <- V(g)$label

    V(g)$x <- sapply( nodedata, function(x) as.numeric( x$kegg.gr.x ) )
    V(g)$y <- -sapply( nodedata, function(x) as.numeric( x$kegg.gr.y ) )

    V(g)$shape <- "rectangle"
    V(g)$shape[comp.nodes] <- "circle"

    V(g)$frame.color <- "black"
    V(g)$frame.color[spare.nodes] <- "gray60"

    V(g)$label.color <- "black"
    V(g)$label.color[spare.nodes] <- "gray60"

    V(g)$label[title.node] <- sub( "TITLE:", "", V(g)$label[title.node] )
    V(g)$label <- sapply( V(g)$label, function(x)
    {
      if( nchar(x) > 25 )
      {
        s <- strsplit(x," ")[[1]]
        return( paste( paste( s[ 1:ceiling(length(s)/2) ], collapse=" "  ),
                       paste( s[ (ceiling(length(s)/2)+1):length(s) ], collapse=" "  ), sep="\n" ) )

      }else return( x )
    })
    V(g)$label[title.node] = paste( V(g)$label[title.node], main, sep="\n")

    V(g)$label.cex <- 0.6
    V(g)$label.cex[title.node] <- 1.2
    V(g)$label.cex[spare.nodes] <- 0.6

    V(g)$color <- "white"
    V(g)$color[title.node] <- "gray75"
    V(g)$color[gene.nodes] <- colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000)[999*(signal.values[gene.nodes]-signal.values.lim[1])/(signal.values.lim[2]-signal.values.lim[1])+1]

    V(g)$size <- 48
    V(g)$size2 <- 20
    V(g)$size[gene.nodes] <- 10
    V(g)$size2[gene.nodes] <- 6
    V(g)$size[spare.nodes] <- 28
    V(g)$size2[spare.nodes] <- 8
    V(g)$size[comp.nodes] <- 8

    E(g)$color <- "gray20"
    E(g)$label.color <- "black"
    E(g)$arrow.size <- 0
    E(g)$arrow.width <- 0
    E(g)$lty <- 3

    E(g)$label <- ""
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "phosphorylation" %in% x ) )] <- "+p"
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "dephosphorylation" %in% x ) )] <- "-p"
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "ubiquitination" %in% x ) )] <- "+u"
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "methylation" %in% x ) )] <- "+m"
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "dissociation" %in% x ) )] <- "#"


    if(!is.null(highlight.genes))
    {
      V(g)$label[which(highlight.genes!="")] <- paste(V(g)[which(highlight.genes!="")]$label," (",highlight.genes[which(highlight.genes!="")],")",sep="")
      V(g)$frame.color[which(highlight.genes!="")] <- "indianred1"
    }
    if(highlight.sinks) V(g)$color[sink.nodes] <- "lightblue2"

    plot(g,asp=0)

    V(g)$label[spare.nodes] <- ""
    V(g)$color[spare.nodes] <- NA
    V(g)$frame.color[spare.nodes] <- NA


    # --> edges
    g2 <- g - E(g)[ which( !sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("activation","expression","reaction") %in% x) ) )  ]
    E(g2)$arrow.size <- 0.4
    E(g2)$arrow.width <- 0.8
    E(g2)$lty <- 1

    par(new=T)
    plot(g2,asp=0)


    # ..> edges
    g2 <- g - E(g)[ which( !sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("indirect effect") %in% x) ) )  ]
    E(g2)$arrow.size <- 0.4
    E(g2)$arrow.width <- 0.8
    E(g2)$lty <- 2

    par(new=T)
    plot(g2,asp=0)


    # --| edges
    g2 <- g - E(g)[ which( !sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("inhibition","repression") %in% x) ) )  ]
    E(g2)$arrow.size <- 0.1
    E(g2)$arrow.width <- 16
    E(g2)$lty <- 1

    par(new=T)
    plot(g2,asp=0)


    # -- edges
    g2 <- g - E(g)[ which( !sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("binding/association","dissociation","compound") %in% x) ) )  ]
    E(g2)$arrow.size <- 0
    E(g2)$arrow.width <- 0
    E(g2)$lty <- 1

    par(new=T)
    plot(g2,asp=0)

  }






  plot.psf.titlepage <- function( psf.object, signal.values )
  {
    layout(matrix(c(1,2,3,4,4,4,5,5,5),3,byrow=TRUE))

    par(mar=c(0,0,0,0))
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
    text(0.05, 0.94, psf.object$attrs$title , cex=2, adj=0)

#     ### Population map of all nodes ###
#     n.map <- matrix(0,preferences$dim.1stLvlSom,preferences$dim.1stLvlSom)
#     pw.genes <- unlist(sapply(psf.object$graph@nodeData@data,function(x)x$genes))
#     pw.metagenes <- som.nodes[names(gene.ids)[which(gene.ids %in% pw.genes)]]
#     n.map[as.numeric(names(table(pw.metagenes)))] <- table(pw.metagenes)
#     n.map[which(n.map==0)] <- NA
#     n.map <- matrix(n.map, preferences$dim.1stLvlSom)
#
#     lim <- c(1,preferences$dim.1stLvlSom) + preferences$dim.1stLvlSom * 0.01 * c(-1, 1)
#     colr <- colramp(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
#                             max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
#                             999 + 1]
#
#     par(mar=c(3,6,3,6))
#     image(matrix(spot.list.overexpression$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE )
#     par(new=TRUE)
#     plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
#          xlab="",ylab="", xaxs="i", yaxs="i", col=colr, main="all genes",
#          cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)
#     title(sub=paste("maximum =", max(n.map,na.rm=TRUE)),line=0)
#     box()

    ### Population map of sink-nodes ###
#     if(length(psf.object$sink.nodes)>0)
#     {
#       n.map <- matrix(0,preferences$dim.1stLvlSom,preferences$dim.1stLvlSom)
#       pw.genes <- unlist(sapply(psf.object$graph@nodeData@data[psf.object$sink.nodes],function(x)x$genes))
#       pw.metagenes <- som.nodes[names(gene.ids)[which(gene.ids %in% pw.genes)]]
#       if(length(pw.metagenes)>0)
#       {
#         n.map[as.numeric(names(table(pw.metagenes)))] <- table(pw.metagenes)
#         n.map[which(n.map==0)] <- NA
#         n.map <- matrix(n.map, preferences$dim.1stLvlSom)
#
#         lim <- c(1,preferences$dim.1stLvlSom) + preferences$dim.1stLvlSom * 0.01 * c(-1, 1)
#         colr <- colramp(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
#                                 max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
#                                 999 + 1]
#
#         par(mar=c(3,6,3,6))
#         image(matrix(spot.list.overexpression$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE )
#         par(new=TRUE)
#         plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
#              xlab="",ylab="", xaxs="i", yaxs="i", col=colr, main="sink node genes",
#              cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)
#         title(sub=paste("maximum =", max(n.map,na.rm=TRUE)),line=0)
#         box()
#
#       } else frame()
#     }


    #       ### Profile of mean sink signals ###
    #     #  hist(  unlist( sapply( psf.results, function(x) sapply(x,function(y) if(length(y$signal.at.sinks)>0) log10(y$signal.at.sinks) else 0 ) ) ),breaks = 40  )
    #       ylim <- c(-5, 5)
    #       par(mar=c(2,7,4,5))
    #
    #       barplot( sapply(signal.values,function(x) mean(log10(x$signal.at.sinks),na.rm=T) ),
    #                 beside=TRUE, col=group.colors, names.arg=rep("",ncol(indata.ensID.m)),
    #                 ylim=ylim, border=if (ncol(indata.ensID.m) < 80) "black" else NA )
    #       mtext(bquote("<log"[10] ~ "s>"), side=2, line=2.5, cex=1.5)
    #
    #       ### Profile of max sink signals ###
    #       ylim <- c(-5, 5)
    #       par(mar=c(6,7,0,5))
    #
    #       bar.coords <- barplot( sapply(signal.values,function(x) max(log10(x$signal.at.sinks),na.rm=T) ),
    #                             beside=TRUE, names.arg=rep("",ncol(indata.ensID.m)),
    #                             col=group.colors, ylim=ylim, border=if (ncol(indata.ensID.m) < 80) "black" else NA )
    #       mtext(bquote("log"[10] ~ "s"[max]), side=2, line=2.5, cex=1.5)
    #
    #       if (ncol(indata.ensID.m)<100)
    #         text(bar.coords, par('usr')[3], labels=colnames(indata.ensID.m), srt=45, adj=c(1.1,1.1), xpd=TRUE)
    #     }
  }



  cat("Writing: PSF signal sheets\n|")
  for( i in 1:50 ) cat(" ");  cat("|\n|");  flush.console()

#   pw=10 #apop
#   pw=14 #bcell rec
#   pw=30 #TCA
  for(pw in 1:length(kegg.collection) )
  {
    filename <- file.path(output.path, paste( make.names(names(kegg.collection)[pw]),".pdf",sep="") )
    pdf(filename, 29.7/2.54, 21/2.54)

    plot.psf.titlepage( psf.object=kegg.collection[[pw]], signal.values=psf.results[[pw]] )

    # node.genes <- lapply(kegg.collection[[pw]]$graph@nodeData@data,function(x)x$genes)
    # node.genes <- sapply(node.genes,head,1)
    # node.genes <- names(gene.ids)[match(node.genes,gene.ids)]

#     node.spots <- rep("",length(node.genes))
#     names(node.spots) <- node.genes
#     for( i in 1:length(spot.list.overexpression$spots) )
#     {
#       node.spots[ which( node.genes %in% spot.list.overexpression$spots[[i]]$genes ) ] <-
#         names(spot.list.overexpression$spots)[i]
#     }
#     names(node.spots) <- names(kegg.collection[[pw]]$graph@nodeData@data)

    par(mfrow=c(1,1),mar=c(0,0,0,0))
    plot.psf.pathway( psf.object = kegg.collection[[pw]])
#     , signal.values = runif(length( psf.results[[pw]]$liver$signal.at.nodes ))/1000,
#                       signal.values.lim = c(-10000,10000), highlight.sinks=TRUE, highlight.genes=node.spots )
# #
# #     for( m in 1:ncol(indata.ensID.m) )
#     {
#       plot.psf.pathway( psf.object = kegg.collection[[pw]], signal.values = log10( psf.results[[pw]][[m]]$signal.at.nodes ),
#                         signal.values.lim = range( log10( sapply( psf.results[[pw]], function(x) x$signal.at.nodes ) ) ),
#                         main = colnames(indata.ensID.m)[m] )
#     }

    dev.off()

    out.intervals <- round( seq( 1, length(kegg.collection), length.out=50+1 ) )[-1]
    cat( paste( rep("#",length( which( out.intervals == pw) ) ), collapse="" ) );	flush.console()
  }
  cat("|\n\n"); flush.console()


}

