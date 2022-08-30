## reading manually selected pathway list
selected_pathway_set <- read.delim(file = "/home/siras/PSFC/inst/extdata/selected_pathways.txt", sep = "\t", stringsAsFactors = F)

## building keggres links data frame for downloading
selected_pathway_set <- data.frame(pathway_code = selected_pathway_set$pathway_code,
                                   kgml_link = paste0("http://rest.kegg.jp/get/", selected_pathway_set$pathway_code, "/kgml"),
                                   png_link = paste0("http://rest.kegg.jp/get/", selected_pathway_set$pathway_code, "/image"),
                                   stringsAsFactors = F
                                   )


## downloading kgmls and png images
setwd("/home/siras/PSFC/inst/extdata/kgmls")
apply(selected_pathway_set, 1, function(x) {
  
  system(paste0("wget ", "-O ", x[1], ".kgml ", x[2]))
  
})
setwd("/home/siras/PSFC/")


setwd("/home/siras/PSFC/inst/extdata/pathway_imgs/")
apply(selected_pathway_set, 1, function(x) {
  
  system(paste0("wget ", "-O ", x[1], ".png ", x[3]))
  
})
setwd("/home/siras/PSFC/")


## parsing kgml files into kegg collecgtion for psf
library(psf)

kegg_collection_new <- generate.kegg.collection.from.kgml(list.files("inst/extdata/kgmls/", full.names = T), sink.nodes = T)

## downlading and parsing group nodes for visualization
### not using anymore ###
## Henry function
library(KEGGgraph)
library(XML)
pathway.kgml = "inst/extdata/kgmls/hsa04151.kgml"

get_group_graphics <- function(pathway.kgml) {
  
  pathway.doc <- try({ xmlTreeParse(pathway.kgml, getDTD = FALSE) }, silent=T )
  if( class(pathway.doc)=="try-error" ) return(NULL)
  
  pathway.doc.r <- xmlRoot(pathway.doc)
  isEntry <- sapply(xmlChildren(pathway.doc.r), xmlName) == "entry"
  
  # attrs <- xmlAttrs(pathway.doc.r)
  # kegg.pathwayinfo <- list(name = attrs[["name"]], title = attrs[["title"]] )
  
  entry=pathway.doc.r[isEntry][[10]]
  kegg.node.info <- lapply( pathway.doc.r[isEntry], function(entry)
  {
    attrs <- xmlAttrs(entry)
    entryID <- getNamedElement(attrs,"id")
    name <- getNamedElement(attrs,"name") # unname(unlist(strsplit(attrs["name"], " ")))
    type <- getNamedElement(attrs,"type")
    
    attrs <- xmlAttrs(xmlChildren(entry)$graphics)
    
    label <- strsplit(getNamedElement(attrs,"name"),", ")[[1]][1]
    label[is.na(label)] = ""
    label <- gsub("[.][.][.]", "", label)
    
    graphics <- list( name = getNamedElement(attrs,"name"), 
                      label = label,
                      x = as.integer(getNamedElement(attrs,"x")), 
                      y = as.integer(getNamedElement(attrs,"y")),  
                      type = getNamedElement(attrs,"type"), 
                      width = as.integer(getNamedElement(attrs,"width")),
                      height = as.integer(getNamedElement(attrs,"height")),        
                      fgcolor = getNamedElement(attrs,"fgcolor"),
                      bgcolor = getNamedElement(attrs,"bgcolor"))
    
    list(entryID = entryID, name = name, type = type, graphics = graphics)
  } )
  names(kegg.node.info) = paste( sapply(kegg.node.info,"[[", "name" ), sapply(kegg.node.info, function(x) x$graphics$x ), sapply(kegg.node.info, function(x) x$graphics$y ) )
  
  
  kegg.node.info <- kegg.node.info[grep("undefined", names(kegg.node.info))]
  
  return(kegg.node.info)
}

group_graphics <- lapply(list.files("inst/extdata/kgmls/", full.names = T), get_group_graphics)

names(group_graphics) <- names(kegg_collection_new)


save(kegg_collection_new, kegg_compounds_to_full_name, entrez_to_symbol, file = "inst/extdata/kegg_collection_new.RData")

## adding group node information to curated pathway collection
for(i in names(edited_pathways_new)){
  edited_pathways_new[[i]]$group_nodes <- kegg_collection_new[[i]]$group_nodes
}

## adding new attributes to old curated pathways
edited_pathways_new <- lapply(edited_pathways_new, function(x) {
  
  graph::edgeDataDefaults(x$graph, attr = "weight") <- 1
  graph::nodeDataDefaults(x$graph, attr = "psf_function") <- "mean"
  
  x
  
})

### removing incorrectly connected map which was clasified as a sink node
edited_pathways_new$MAPK_signaling_pathway$graph <- graph::removeEdge(from = "47", to = "42", graph = edited_pathways_new$MAPK_signaling_pathway$graph)
edited_pathways_new$MAPK_signaling_pathway$graph <- psf::set.edge.impacts(edited_pathways_new$MAPK_signaling_pathway$graph)
edited_pathways_new$MAPK_signaling_pathway$order <- psf::order.nodes(edited_pathways_new$MAPK_signaling_pathway$graph)
edited_pathways_new$MAPK_signaling_pathway$sink.nodes <- psf::determine.sink.nodes(edited_pathways_new$MAPK_signaling_pathway)



save(edited_pathways_new, kegg_compounds_to_full_name, entrez_to_symbol, file = "inst/extdata/edited_pathways_new.RData")

