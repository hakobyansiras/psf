# load("../oposSOM/Tissues.RData")
# load("../psf_scripts/KEGG.Collection.Rdata")
# # kgml = "data_to_be_deleted/KGML/hsa04662.xml"
# # g = parse.KGML(kgml)
# # indata = env$indata
# # indata.FC = 10^indata
# # entrez.FC = ensembl.id.conversion(indata.FC)
# # gg = map.gene.data(g,entrez.FC)
# 
# compute.flows.from.env <- function(env, g, node.ordering){
#   indata = env$indata
#   indata.FC = 10^indata
#   entrez.FC = ensembl.id.conversion(indata.FC)
#   g = map.gene.data(g, entrez.FC)
#   psf.results = psf.flow(g, node.ordering)
# }
# 
# compute.flows.for.KEGG.Collection <- function(env){
#   if(is.null(env))
#     stop("null env argument")
# 
#   cat("PSF computation started.\n")
#   # cat("Loading KEGG Collection.\n")
# 
#   # load("KEGG.Collection.Rdata") # TBA: provide the location of the collection
# 
# 
# #   if(!("kegg.collection" %in% ls()))
# #     stop("Problem loading KEGG Collection")
#   if(length(kegg.collection) == 0)
#     stop("Empty KEGG Collection")
#   cat("KEGG.Collection with ", length(kegg.collection), " pathways successfully loaded\n")
# 
#   #genewise normalized data: considered as logFC.
#   indata = env$indata
# 
#   #anti-log the indata to get the FC values
#   indata.FC = 10^indata
#   cat("Converting ENSEMBL gene ids to Entrez gene ids: this may take a while.\n")
#   entrez.FC <<- ensembl.id.conversion(indata.FC)
# 
#   #Starting the psf calculation.
#   # create psf.results list
#   #
#   # for each pathway in the collection {
#   #   map the gene expression values onto the pathway
#   #   compute psf values for that pathway
#   #   add to psf.results
#   # }
# 
# 
#   # genes are kept in rownames
#   # gene ID conversion
#   #
#   cat("PSF calculations on the collection")
#   psf.results.collection = list()
#   for(i in 1:length(kegg.collection)){
#     pathway = names(kegg.collection)[i]
#     if(!length(kegg.collection[[i]])==0){
#       cat("Current pathway: ", pathway, "\n")
#       g = map.gene.data(kegg.collection[[i]]$graph, entrez.FC)
#       psf.results.collection$pathway = psf.flow(g, kegg.collection[[i]]$order)
#     }
#   }
#   psf.results.collection[[1]]$signal.at.sink
# 
# 
# }
