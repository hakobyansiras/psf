#' Download a KEGG pathway from the web
#' @param id KEGG identifier of the pathway
#' @param dir The directory where the downloaded file should be saved
#' @return The filepath of the downloaded pathway or NULL the pathway was not downloaded
#'
download.KGML <- function(id, dir){
  if (!dir.exists(dir)){
    if(!dir.create(dir))
      stop("Directory ", dir, " could not be created.")
  }

  #! replace with tryCatch
  file.content <- try(KEGGREST::keggGet(id, option = c("kgml")), silent = T)
  if (class(file.content)!="try-error") {
    doc = XML::xmlParseDoc(file.content, asText = T)
    kgml.file = file.path(dir,paste(id, ".xml", sep=""))
    XML::saveXML(doc,file=kgml.file)
    return(kgml.file)
  } else {
    stop("The KGML file ", id, " could not be downloaded")
  }
  return(file)
}

#' Parse KGML to a graph object
#'
#' @param kgml The KGML file path
#' @return GraphNEL object containing the parsed pathway
#'
parse.KGML <- function(kgml){
  if(!file.exists(kgml))
    stop("Provided KGML file: ", kgml, ", does not exist")

  doc <-XML::xmlParseDoc(kgml)
  r <- XML::xmlRoot(doc)
  root.children = XML::xmlChildren(r)


  # Reading the KEGG attributes of <entry> nodes and storing them in nodelist.mat.

  entry.nodes = root.children[which(names(root.children) == "entry")]

  entry.attr.names = c("kegg.id","kegg.name", "kegg.type", "kegg.link")
  names(entry.attr.names) = c("id", "name", "type", "link")
  graphics.attr.names = c("kegg.gr.name","kegg.gr.fgcolor","kegg.gr.bgcolor","kegg.gr.type", "kegg.gr.x", "kegg.gr.y", "kegg.gr.width", "kegg.gr.height")
  names(graphics.attr.names) = c("name", "fgcolor", "bgcolor", "type", "x", "y", "width", "height")

  nodelist.mat = matrix("NA", nrow=length(entry.nodes), ncol = length(c(entry.attr.names, graphics.attr.names)) + 1)
  colnames(nodelist.mat) = c(entry.attr.names, graphics.attr.names, "components")


  i = 1
  for (entry in entry.nodes){
    entry.attrs = XML::xmlAttrs(entry)
    for(name in names(entry.attr.names)){
      if(name %in% names(entry.attrs))
        nodelist.mat[i, entry.attr.names[[name]]] = as.character(entry.attrs[[name]])
    }
    entry.children = XML::xmlChildren(entry)
    graphics = entry.children$graphics
    graphics.attrs = XML::xmlAttrs(graphics)
    for(name in names(graphics.attr.names)){
      if(name %in% names(graphics.attrs))
        nodelist.mat[i, graphics.attr.names[[name]]] = as.character(graphics.attrs[[name]])
    }

    if(entry.attrs[["type"]] == "group"){
      components = entry.children[which(names(entry.children) == "component")]
      if(length(components) > 0){
        componentids = ""
        for(component in components){
          id = as.character(XML::xmlGetAttr(component, "id"))
          componentids = paste(componentids, id, sep=";")
        }
      }
      nodelist.mat[i, "components"] = componentids
    }

    if(graphics.attrs[["type"]] == "line"){
      coords = XML::xmlGetAttr(graphics,"coords")
      # cat("parsing ", coords)
      #assume coords has the form [x1, y1, x2, y2]
      coords = as.numeric(unlist(strsplit(coords, split =",")))
      xmean = base::mean(c(coords[1], coords[3]))
      ymean = base::mean(c(coords[2], coords[4]))
      nodelist.mat[i,"kegg.gr.x"] = xmean
      nodelist.mat[i,"kegg.gr.y"] = ymean
      # cat(" xmean=", xmean, " ymean=",ymean, "\n")
    }

    i = i + 1
  }
  rownames(nodelist.mat) = nodelist.mat[,"kegg.id"]

  # Building a graphNEL object with the parsed entry nodes.

  g = graph::graphNEL(nodes = c(rownames(nodelist.mat)), edgemode = "directed")

  graph::nodeDataDefaults(g,attr="genes") <- list()
  graph::nodeDataDefaults(g,attr="expression") <- 1
  graph::nodeDataDefaults(g,attr="signal") <- 1
  graph::nodeDataDefaults(g,attr="type") <- "gene"
  graph::nodeDataDefaults(g,attr="label") <- "NA"
  graph::nodeDataDefaults(g,attr="components") <- "NA"
  for(attr in entry.attr.names){
    graph::nodeDataDefaults(g, attr = attr) <- NA
  }
  for(attr in graphics.attr.names){
    graph::nodeDataDefaults(g, attr = attr) <- NA
  }

  for(attr in entry.attr.names){
    graph::nodeData(g, graph::nodes(g), attr = attr) <- nodelist.mat[,attr]
  }

  for(attr in graphics.attr.names){
    graph::nodeData(g, graph::nodes(g), attr = attr) <- nodelist.mat[,attr]
  }
  graph::nodeData(g, graph::nodes(g), "components") = nodelist.mat[,"components"]
  graph::nodeData(g, graph::nodes(g), "type") = nodelist.mat[,"kegg.type"]

  # The component gene entrez ids are kept in the kegg.name attribute,
  # with "hsa:" prefix. We split this attribute value, and keep the
  # component gene entrez ids in the "genes" attribute for each node.

  # The kegg.gr.name contains component symbols, we take the first symbol
  # as the node label (used for graph visualization)



  for(node in graph::nodes(g)){
    type = nodelist.mat[node, "kegg.type"]
    if(type == "gene") {
      kegg.name = nodelist.mat[node, "kegg.name"]
      if(!is.na(kegg.name)){
        entrezIds = unlist(base::strsplit(kegg.name, "hsa:", fixed = T))
        entrezIds = unlist(base::strsplit(entrezIds, " ", fixed = T))
        na.ind = which(entrezIds == "")
        if(length(na.ind) > 0)
          entrezIds = entrezIds[-na.ind]
      }
      if(length(entrezIds) > 0 )
        graph::nodeData(g, node, "genes") = list(entrezIds)
    }

    kegg.gr.name = nodelist.mat[node, "kegg.gr.name"]
    gr.names = unlist(strsplit(kegg.gr.name, ", "))
    na.ind = which(gr.names == "")
    if(length(na.ind) > 0)
      gr.names = gr.names[-na.ind]
    label = gr.names[1]
    graph::nodeData(g, node, "label") = label

  }

  # Adding edges to the graph

  edge.attrs = list(impact="impact", type="type", subtype1="subtype1",
                    subtypeValue1 = "subtypeValue1",
                    subtypeValue2 = "subtypeValue2", subtype2 = "subtype2")
  graph::edgeDataDefaults(g,attr=edge.attrs$impact) <- 1
  graph::edgeDataDefaults(g,attr=edge.attrs$type) <- NA
  graph::edgeDataDefaults(g,attr=edge.attrs$subtype1) <- NA
  graph::edgeDataDefaults(g,attr=edge.attrs$subtypeValue1) <- NA
  graph::edgeDataDefaults(g,attr=edge.attrs$subtype2) <- NA
  graph::edgeDataDefaults(g,attr=edge.attrs$subtypeValue2) <- NA


  relations = root.children[which(names(root.children) == "relation")]
  for(relation in relations){
    entry1 = as.character(XML::xmlGetAttr(relation, "entry1"))
    entry2 = as.character(XML::xmlGetAttr(relation, "entry2"))
    type = as.character(XML::xmlGetAttr(relation, "type"))
    g = psf::add.kegg.edge(entry1, entry2, type, subtype="n/a", edge.attrs, g)
    relation.children = XML::xmlChildren(relation)
    subtypes = relation.children[which(names(relation.children) == "subtype")]
    if(length(subtypes) >= 1){
      subtype1atts = XML::xmlAttrs(subtypes[[1]])
      if(entry2 %in% graph::edges(g)[[entry1]]){
        graph::edgeData(g, entry1, entry2, edge.attrs$subtype1) = as.character(subtype1atts[["name"]])
        graph::edgeData(g, entry1, entry2, edge.attrs$subtypeValue1) = as.character(subtype1atts[["value"]])
        if(length(subtypes) == 2){
          subtype2atts = XML::xmlAttrs(subtypes[[2]])
          graph::edgeData(g, entry1, entry2, edge.attrs$subtype2) = as.character(subtype2atts[["name"]])
          graph::edgeData(g, entry1, entry2, edge.attrs$subtypeValue2) = as.character(subtype2atts[["value"]])
        }
      } else {
        cat("Edge ", entry1, ":", entry2, " not found\n")
      }
    } else {
      cat("No subtypes found for edge ", entry1, "|", entry2 , "\n")
    }
  }

  #####################################################
  ###            reactions                      #######
  #####################################################
  reactions = root.children[which(names(root.children) == "reaction")]
  for(reaction in reactions){
    reaction.id = XML::xmlGetAttr(reaction, "id")
    reaction.name = XML::xmlGetAttr(reaction, "name")
    reaction.entries = XML::xmlChildren(reaction)
    subst.ind = which(names(reaction.entries) == "substrate")
    prod.ind = which(names(reaction.entries) == "product")
    ### edge for substrate - enzyme
    if(length(subst.ind) > 1)
      cat("\nreaction ", reaction.name, " had " , length(subst.ind), " substrates\n")
    for(s.id in subst.ind){
      substrate = reaction.entries[[s.id]]
      substrate.id = XML::xmlGetAttr(substrate, "id")
      entry1 = substrate.id
      entry2 = reaction.id
      g = psf::add.kegg.edge(entry1, entry2, "ECrel", "reaction", edge.attrs, g)
    }

    ### edge for enzyme - product
    if(length(prod.ind) > 1)
      cat("\nreaction ", reaction.name, " had " , length(prod.ind), " products\n")
    for(p.id in prod.ind){
      product = reaction.entries[[p.id]]
      product.id = XML::xmlGetAttr(product, "id")
      entry1 = reaction.id
      entry2 = product.id
      g = psf::add.kegg.edge(entry1, entry2, "ECrel", "reaction", edge.attrs, g)
    }
  }


  #####################################################
  ###            end of reactions               #######
  #####################################################

  ## handle group nodes: connect component nodes to each other
  ## and to the rest of the graph and remove group nodes
  for(node in graph::nodes(g)){
    if(graph::nodeData(g, node, "type") == "group"){
      components = as.character(graph::nodeData(g, node, "components"))
      components = unlist(strsplit(components, split = ";"))
      na.ind = which(components == "")
      if(length(na.ind) > 0)
        components = components[-na.ind]
      g = process.groupNode(node,components, g, edge.attrs)
      g = graph::removeNode(node, g)
      # cat("Group node processed: ", node, "\n")
    }
  }

  ## turn gene->gene interactions to gene->compound->gene interactions where applicable
  g = process.compounds(g, edge.attrs)

  ## check for binding interaction directions -> if guessed to be wrong, reverse them
  g = correctEdgeDirections(g, edge.attrs)

  ## set edge impacts
  g = set.edge.impacts(g)


  return(g)
}

add.kegg.edge.mut <- function(entry1, entry2, subtype, g){
  edge.attrs = list(impact="impact", type="type", subtype1="subtype1",
                    subtypeValue1 = "subtypeValue1",
                    subtypeValue2 = "subtypeValue2", subtype2 = "subtype2")
  g = add.kegg.edge(entry1, entry2, "PPrel", subtype, edge.attrs, g)

  if(subtype == "inhibition")
    graph::edgeData(g, entry1, entry2, attr = "impact") <- -1
  else
    graph::edgeData(g, entry1, entry2, attr = "impact") <- 1
  return(g)
}

add.kegg.edge <- function(entry1, entry2, type, subtype, edge.attrs, g){
  if(is.null(type) || is.na(type))
    stop("null type for edge ", entry1, "|", entry2, "\n")
  if(is.null(subtype) || is.na(subtype))
    stop("null subtype for edge ", entry1, "|", entry2, "\n")
  if(!(entry1 %in% graph::nodes(g)))
    stop("entry1 " + entry1 + " not in graph")
  if(!(entry2 %in% graph::nodes(g)))
    stop("entry2 " + entry2 + " not in graph")
  names(type) = NULL
  names(subtype) = NULL
    if(!is.null(entry1) & !is.null(entry2)){
    if(graph::isAdjacent(g, entry1, entry2)){
      cat("Edge ", entry1, "|", entry2, "already exists\n")
      if(is.na(graph::edgeData(g, entry1, entry2, edge.attrs$type)))
        graph::edgeData(g, entry1, entry2, edge.attrs$type) = type
      if(is.na(graph::edgeData(g, entry1, entry2, edge.attrs$subtype1)))
        graph::edgeData(g, entry1, entry2, edge.attrs$subtype1) = subtype
    }
    else{
      g = graph::addEdge(entry1, entry2, g)
      graph::edgeData(g, entry1, entry2, edge.attrs$type) = type
      graph::edgeData(g, entry1, entry2, edge.attrs$subtype1) = subtype
    }
  }
  return(g)
}

set.edge.impacts <- function(g, subtype1.attr = "subtype1",
                             subtype2.attr = "subtype2", impact.attr.name = "impact"){
  inhibition.ind.1 <-grep("inhibition",  unlist(graph::edgeData(g, attr = subtype1.attr)))
  inhibition.ind.2 <-grep("inhibition",  unlist(graph::edgeData(g, attr = subtype2.attr)))
  inhibition.ind = union(inhibition.ind.1, inhibition.ind.2)
  repression.ind.1 <-grep("repression",  unlist(graph::edgeData(g, attr = subtype1.attr)))
  repression.ind.2 <-grep("repression",  unlist(graph::edgeData(g, attr = subtype2.attr)))
  repression.ind = union(repression.ind.1, inhibition.ind.2)

  inhibition.ind <- union(inhibition.ind,repression.ind)
  if (length(inhibition.ind)>0){
    from.to = strsplit(names(unlist(graph::edgeData(g,attr=impact.attr.name))[inhibition.ind]),
                       "|", fixed=T)
    for (i in 1:length(from.to)){
      tmp <- from.to[[i]]
      graph::edgeData(g, tmp[1], tmp[2], attr = impact.attr.name) <- -1
    }
  }

  return(g)
}


get.pathway.attrs <- function(kgml){
  if(!file.exists(kgml))
    stop("Provided KGML file: ", kgml, ", does not exist")

  doc <-XML::xmlParseDoc(kgml)
  r <- XML::xmlRoot(doc)

  #Reading pathway attributes

  pathway.name = XML::xmlGetAttr(r,"name")
  pathway.title = XML::xmlGetAttr(r, "title")
  pathway.link =  XML::xmlGetAttr(r, "link")
  pathway.image =  XML::xmlGetAttr(r, "image")

  pathway.attrs = list(name = pathway.name, title=pathway.title, image = pathway.image, link = pathway.link)
  return(pathway.attrs)
}

#' Extend the group node to its component gene nodes
#'
#' @param groupNode The group node to be processed
#' @param components The vector of its component gene nodes
#' @param g The pathway graph of type GraphNEL
#' @param edge.attrs The list of edge attr
#'
#' @details The all the edges to and from the group node are
process.groupNode <- function(groupNode, components, g, edge.attrs){

  inComponents = vector(mode = "character") #components receiving incoming connections
  outComponents = vector(mode = "character")  #components giving outgoing connections
  biComponents = vector(mode = "character") #components in the middle
  incomingNodes = vector(mode = "character") #non-group nodes with incoming connections to group components
  outgoingNodes = vector(mode = "character")  #non-group nodes with receiving connections from group components
  incomingGroupNodes = vector(mode = "character") #non-group nodes with incoming connections to the group node itself
  outgoingGroupNodes = vector(mode = "character") #non-group nodes recieving connections from the group node itself


  for(node in graph::nodes(g)){
    if (groupNode %in% graph::edges(g)[[node]]) {
      if(!(node %in% incomingGroupNodes))
        incomingGroupNodes = c(incomingGroupNodes,node)
    } else if(node %in% graph::edges(g)[[groupNode]]) {
      if (!(node %in% outgoingGroupNodes))
        outgoingGroupNodes = c(outgoingGroupNodes, node);
    }
    size = length(components)
    for (component in components) {
      if (component %in% graph::edges(g)[[node]]) {
        if (!(node %in% incomingNodes))
          incomingNodes = c(incomingNodes, node)
        if (component %in% outComponents) {
          if (!(component %in% biComponents))
            biComponents = c(biComponents, component)
          outComponents = outComponents[-which(outComponents == component)]
        } else if (!(component %in% inComponents))
          inComponents = c(inComponents, component)
      } else if (node %in% graph::edges(g)[[component]]) {
        if (component %in% inComponents) {
          if (!(component %in% biComponents))
            biComponents = c(biComponents, component)
          inComponents = inComponents[-which(inComponents == component)]
        } else if (!(component %in% outComponents))
          outComponents = c(outComponents, component)
        if (!(node %in% outgoingNodes))
          outgoingNodes = c(outgoingNodes, node)
      }
    }
  }

  for (component in inComponents)
    components = components[-which(components == component)]
  for (component in outComponents)
    components = components[-which(components == component)]
  for (component in biComponents)
    components = components[-which(components == component)]

  prevNode = NULL

  if(length(biComponents) > 0) {
    if(length(inComponents) > 0) {
      if(length(outComponents)==0){
        for(biNode in biComponents) {
          if(!(biNode %in% outComponents))
            outComponents = c(outComponents, biNode)
        }
      }
    } else for(biNode in biComponents) {
      if(!(biNode %in% inComponents))
        inComponents = c(inComponents, biNode)
    }
  }

  #   Make the first node the innode. Remove all other incomming edges and move them to
  #   the innode.

  if(length(inComponents) > 0) {
    inNode = inComponents[1]
    prevNode = inNode
    inComponents = inComponents[-1]
    if(length(inComponents) > 0)
      for(componentNode in inComponents) {
        for(incomingNode in incomingNodes)
          if(componentNode %in% graph::edges(g)[[incomingNode]]) {
            g = redirectEdge(incomingNode, componentNode, incomingNode, inNode, g, edge.attrs)
          }
        g = psf::add.kegg.edge(prevNode, componentNode,type="ECrel", subtype="compound", edge.attrs, g)
        prevNode = componentNode
      }

  } else if(length(components) > 0) {
    inNode = components[1]
    components = components[-1]
    prevNode = inNode;
  } else if(length(outComponents) > 0) {
    inNode = outComponents[1]
  }


  while(length(components) > 0) {
    if(!is.null(prevNode))
      g = psf::add.kegg.edge(prevNode, components[1],"ECrel", "compound", edge.attrs, g)
    prevNode = components[1]
    components = components[-1]
  }
  #   Make the first node the outNode. Remove all the edges from other
  #   components in outComponents and redirect them from the outNode.

  if(length(outComponents) > 0) {
    outNode = outComponents[1]
    outComponents = outComponents[-1]
    if(length(outComponents) > 0)
      for(componentNode in outComponents) {
        for(outgoingNode in outgoingNodes)
          if(outgoingNode %in% graph::edges(g)[[componentNode]]) {
            g = redirectEdge(componentNode, outgoingNode, outNode, outgoingNode, g, edge.attrs)
          }
        if (!is.null(prevNode))
          g = psf::add.kegg.edge(prevNode, componentNode, "ECrel","compound", edge.attrs, g)
        prevNode = componentNode
      }
    g = psf::add.kegg.edge(prevNode, outNode, "ECrel","compound", edge.attrs, g)
  } else
    outNode = prevNode;

  #   Add the edges to the group node to the inNode and outNode.
  #   Remove group node.

  if(length(incomingGroupNodes) != 0)
    for(incomingGroupNode in incomingGroupNodes) {
      g = redirectEdge(incomingGroupNode, groupNode, incomingGroupNode, inNode, g, edge.attrs)
    }
  if(length(outgoingGroupNodes) != 0)
    for(outgoingGroupNode in outgoingGroupNodes) {
      g = redirectEdge(groupNode, outgoingGroupNode, outNode, outgoingGroupNode, g, edge.attrs)
    }
  for(outGoingNode in outgoingNodes) {
    if(outGoingNode %in% graph::edges(g)[[inNode]]) {
      g = redirectEdge(inNode, outGoingNode, outNode, outGoingNode, g, edge.attrs)
    }
  }

  for(inComingNode in incomingNodes){
    if(outNode %in% graph::edges(g)[[inComingNode]]) {
      g = redirectEdge(inComingNode, outNode, inComingNode, inNode, g, edge.attrs)
    }
  }
  return(g)
}

#' Add an edge between two nodes if no such edge exists,
#' and if no reverse edge exists
#'
#' @param sourceNode The source node
#' @param targetNode The target node
#' @param g The graph of class graphNEL
#'
#' @return g The modified graph with added (or not) edge
# addEdgeSafe <- function(sourceNode, targetNode, g) {
#   if(!(targetNode %in% graph::edges(g)[[sourceNode]]))
#     # if (!(sourceNode %in% graph::edges(g)[[targetNode]])) {
#     g = graph::addEdge(sourceNode, targetNode, g)
#   # }
#   return(g)
# }


#' Redirect the edge to new source and target nodes
#'
#' @param prevSource The original source node
#' @param prevTarget The original target node
#' @param newSource The new source node
#' @param newTarget The new target node
#' @param g The pathway graph of graphNEL class
#' @param edge.attrs The list of edge attributes in the graph
#'
#' @details A new edge between the new source and target nodes will be added to the graph, the previous edge will be removed, and its attributes will be transfered to the new edge.
redirectEdge <- function(prevSource, prevTarget,newSource, newTarget, g, edge.attrs) {
  if(prevTarget %in% graph::edges(g)[[prevSource]]) {
    prev.edge.attrs = list()
    for(attr in edge.attrs){
      prev.edge.attrs[[attr]] = graph::edgeData(g, prevSource, prevTarget, attr)
    }

    g = graph::removeEdge(prevSource, prevTarget, g)
    g = psf::add.kegg.edge(newSource, newTarget, prev.edge.attrs$type, prev.edge.attrs$subtype1, edge.attrs, g)
    if(newTarget %in% graph::edges(g)[[newSource]])
      for(attr in edge.attrs){
        graph::edgeData(g, newSource, newTarget, attr) = prev.edge.attrs[[attr]]
      }

  }
  return(g)
}

change.edge.type <- function(from, to, subtype, g) {
  g = psf::remove.edge(from, to, g)
  g = psf::add.kegg.edge.mut(from, to, subtype, g)
  return(g)
}

remove.edge <- function(from, to, g){
  g = graph::removeEdge(from, to, g)
  return(g)
}

#' Process protein compound interactions
#'
#'@param g The pathway graph of graphNEL class
#'@param edge.attrs The list of edge attributes in the pathway graph
#'
#'@details This function will turn gene -> gene interactions that occur via a compound node into gene -> compound -> gene interactions.
process.compounds <- function(g, edge.attrs) {

  compoundRelations = vector(mode="character") #Relations to be removed
  newRelations = vector(mode="character") #Relations to be added
  for (snode in names(graph::edges(g))) {
    tnodes = graph::edges(g)[[snode]]

    if(length(tnodes) != 0){
      for(tnode in tnodes){
        if (graph::edgeData(g, snode, tnode, edge.attrs$subtype1) == "compound" ||
            graph::edgeData(g, snode, tnode, edge.attrs$subtype2) == "compound") {
          # cat("Compound processing for edge: ", snode, tnode, "\n")
          prev.edge.attrs = list()
          for(attr in edge.attrs){
            prev.edge.attrs$attr = graph::edgeData(g, snode, tnode, attr)
          }


          compound = NULL
          if(graph::edgeData(g, snode, tnode, edge.attrs$subtype1) == "compound")
            compound = graph::edgeData(g, snode, tnode, edge.attrs$subtypeValue1)
          else
            compound = graph::edgeData(g, snode, tnode, edge.attrs$subtypeValue2)
          compound = as.character(compound)

          if(!is.null(compound) && compound %in% graph::nodes(g)) {
            g = psf::add.kegg.edge(snode, compound, "ECrel", "compound", edge.attrs, g)
            # cat("Compound processed: ", snode, ":", compound, " edge added", "\n")
            g = psf::add.kegg.edge(compound, tnode, "ECrel", "compound", edge.attrs, g)
            # cat("Compound processed: ", compound, ":", tnode, " edge added", "\n")

            if(graph::edgeData(g, snode, tnode, edge.attrs$subtype1) == "compound"){
              otherSubtype = graph::edgeData(g, snode, tnode, edge.attrs$subtype2)
              otherSubtypeValue = graph::edgeData(g, snode, tnode, edge.attrs$subtypeValue2)
            } else{
              otherSubtype = graph::edgeData(g, snode, tnode, edge.attrs$subtype1)
              otherSubtypeValue = graph::edgeData(g, snode, tnode, edge.attrs$subtypeValue1)
            }


            graph::edgeData(g, snode, compound, edge.attrs$subtype2) = otherSubtype
            graph::edgeData(g, snode, compound, edge.attrs$subtypeValue2) = otherSubtypeValue

            g = graph::removeEdge(snode, tnode, g)
            # cat("Compound processed: ", snode, ":", tnode, " edge removed", "\n")

          }
        }
      }
    }
  }
  return(g)
}

#' Guess wrong directed binding interactions and reverse them
#' @param g The graph object of graphNEL class
#' @param edge.attrs The list of edge attributes in the graph
#'
#' @return g The modified graph
correctEdgeDirections <- function(g, edge.attrs) {
  #If relations are from higher Id to lower, leave it as it is.
  # Otherwise, check, if the second node lies lefter of higher than the first one, reverse the edge.

  for(snode in names(graph::edges(g))) {
    for(tnode in graph::edges(g)[[snode]]){
      if (graph::edgeData(g, snode, tnode, edge.attrs$subtype1) == "binding" ||
          graph::edgeData(g, snode, tnode, edge.attrs$subtype2) == "binding")
        if (as.numeric(snode) < as.numeric(tnode)){
          if (isReverseDirection(snode, tnode, g)) {
            g = reverseEdges.add(snode, tnode, g, edge.attrs)
            # cat("Binding interaction directions: Edge ", snode, ":", tnode, " was reversed.\n")
          }
        }
    }
  }

  return(g)

}

isReverseDirection <- function(snode, tnode, g) {
  sx = graph::nodeData(g, snode, "kegg.x")
  sy = graph::nodeData(g, snode, "kegg.y")
  tx = graph::nodeData(g, tnode, "kegg.x")
  ty = graph::nodeData(g, snode, "kegg.y")
  if (sx - tx != 0)
    return(sx > tx)
  else
    return(sy > ty)
}

reverseEdge <- function(snode, tnode, g, edge.attrs) {
  return(redirectEdge(snode, tnode, tnode, snode, g, edge.attrs))
}

remove.disconnected.nodes <- function(g){
  edges = graph::edges(g)
  connected.nodes = vector()
  for(snode in names(edges)){
    if(length(edges[[snode]]) > 0){
      if(!(snode %in% connected.nodes))
        connected.nodes = c(connected.nodes, snode)
      for(tnode in edges[[snode]]){
        if(!(tnode %in% connected.nodes))
          connected.nodes = c(connected.nodes, tnode)
      }
    }
  }
  all.nodes = graph::nodes(g)
  disconnected.nodes = setdiff(all.nodes, connected.nodes)
  for(node in disconnected.nodes){
    g = graph::removeNode(node, g)
  }
  return(g)
}

#' Plots the pathway with KEGG layout based on x and y coordinates taken from KGML file.
plot.pathway <- function(g, sink.nodes = NULL){
  igr = igraph::igraph.from.graphNEL(g)
  node.labels = graph::nodeData(g, graph::nodes(g), "label")
  igraph::V(igr)$name = node.labels
  #for now
  # igraph::V(igr)$label = graph::nodeData(g, graph::nodes(g), "kegg.id")
  igraph::V(igr)$ID = graph::nodes(g)

  coords = matrix(data = NA, nrow = length(graph::nodes(g)), ncol = 2)
  coords.x = as.numeric(unlist(graph::nodeData(g, graph::nodes(g), "kegg.gr.x")))
  coords.y = as.numeric(unlist(graph::nodeData(g, graph::nodes(g), "kegg.gr.y")))
  coords.y = (max(coords.y) - coords.y)
  igraph::V(igr)$x = coords.x
  igraph::V(igr)$y = coords.y

  node.shapes = rep_len("circle", length(graph::nodes(g)))
  names(node.shapes) = graph::nodes(g)
  for(node in graph::nodes(g)){
    if(graph::nodeData(g, node, "kegg.type") == "gene")
      node.shapes[node] = "rectangle"
  }

  lapply(names(g@edgeData@data), function(x) {
    get.edge.type(g,
                  strsplit(x, split="|", fixed = T)[[1]][1],
                  strsplit(x, split="|", fixed = T)[[1]][2])
  })


  igraph::V(igr)$vertex.shape = node.shapes
  if(!is.null(sink.nodes)){
    igraph::V(igr)$color <- ifelse((igraph::V(igr)$ID %in% sink.nodes), "red", "blue")
  }

  igraph::tkplot(igr, vertex.shape = igraph::V(igr)$vertex.shape)
  plot(igr, vertex.shape = igraph::V(igr)$shape)

}

generate.kegg.collection <- function(pathway.id.list, out.dir, sink.nodes = T){
  kegg.collection = list()
  for(id in pathway.id.list){
    kgml = download.KGML(id, out.dir)
    cat("downloaded: ", id, "\n")

    if(!is.null(kgml) && file.exists(kgml)){
      pathway.attrs = get.pathway.attrs(kgml)
      g = parse.KGML(kgml)
      order = order.nodes(g)
      title = pathway.attrs$title
      kegg.collection[[title]] = list(graph=g, order=order, attrs=pathway.attrs)
    }
    cat("parsed: " , kgml, "\n")
  }
  if(sink.nodes){
    kegg.collection = determine.sink.nodes.for.collection(kegg.collection)
  }
  return(kegg.collection)
}

generate.kegg.collection.from.kgml <- function(kgml.files, out.dir, sink.nodes = T){
  kegg.collection = list()
  for(kgml in kgml.files){
    cat("kgml: ", kgml, "\n")

    if(!is.null(kgml) && file.exists(kgml)){
      pathway.attrs = get.pathway.attrs(kgml)
      g = parse.KGML(kgml)
      order = order.nodes(g)
      title = pathway.attrs$title
      kegg.collection[[title]] = list(graph=g, order=order, attrs=pathway.attrs)
    }
    cat("parsed: " , kgml, "\n")
  }
  if(sink.nodes){
    kegg.collection = determine.sink.nodes.for.collection(kegg.collection)
  }
  return(kegg.collection)
}

convert.ids.in.kegg.collection <- function(kegg.collection){
  ### ensembl.id.conversion ###

  all.entrez.ids <- unique( unlist( sapply( kegg.collection, function(x)
  {
    d <- x$graph@nodeData@data
    sapply(d,function(y) y$genes )
  }) ) )

  mart <- biomaRt::useMart("ensembl")
  mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = mart)
  biomart.table <- biomaRt::getBM(c("entrezgene","ensembl_gene_id", "hgnc_symbol"),
                                  values = all.entrez.ids, mart = mart, filters = "entrezgene",
                                  uniqueRows = T)


  entrez.to.ensembl <- biomart.table[ match( unique(biomart.table[,1]), biomart.table[,1] ), 2 ]
  entrez.to.symbol <- biomart.table[ match( unique(biomart.table[,1]), biomart.table[,1] ), 3 ]
  names(entrez.to.ensembl) = unique(biomart.table[,1])
  names(entrez.to.symbol) = unique(biomart.table[,1])


  kegg.collection <- lapply( kegg.collection, function(x)
  {
    if( length(x) > 0 )
    {
      graph::nodeDataDefaults(x$graph,attr="entrez.genes") <- list()
      graph::nodeDataDefaults(x$graph,attr="gene.symbols") <- list()
      x$graph@nodeData@data <- lapply( x$graph@nodeData@data, function(y)
      {
        y$entrez.genes <- y$genes
        y$genes <- na.omit(entrez.to.ensembl[y$entrez.genes])
        y$gene.symbols <- na.omit(entrez.to.symbol[y$entrez.genes])
        return(y)
      })
    }
    return(x)
  })
  return(kegg.collection)
}

order.nodes <- function(g){
  g = Rgraphviz::layoutGraph(g)
  nodeY <-graph::nodeRenderInfo(g)$nodeY
  node.rank <-base:::rank(nodeY, ties.method="min")
  node.order <- base:::sort(node.rank,decreasing=T)
  order = list("node.order" = node.order, "node.rank" = node.rank)
  return(order)
}

