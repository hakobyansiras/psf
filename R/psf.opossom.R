perform.psf <- function(env, entrez.fc = NULL,kegg.collection=NULL, calculate.significance = T, bst.steps = 200, split = TRUE){
  cat("Peforming PSF calculations.\n")
  if(is.null(env))
    stop("null env argument")
  if(is.null(kegg.collection)){
    cat("Loading kegg collection\n")
    kegg.collection.dir = "inst/extdata/kegg.collection.RData"
    cat("Trying ", kegg.collection.dir, "\n")
    load(kegg.collection.dir)
    if(is.null(kegg.collection))
      stop("kegg collection not found in ", kegg.collection.dir, "\n")
  }
  if(length(kegg.collection) == 0)
    stop("Empty kegg.collection")
  cat("kegg.collection with ", length(kegg.collection), " pathways successfully loaded\n")

  if(is.null(entrez.fc)){
    cat("Extracting data\n")
    #genewise normalized data: considered as logFC.
    indata = env$indata

    #anti-log the indata to get the FC values
    indata.fc = 10^indata

    cat("Converting ENSEMBL IDs to Entrez gene ids\n")
    entrez.fc = ensembl.id.conversion(indata.fc)
  }
  cat("\nDatasets loaded. Performing PSF calculations\n\n")
  psf.results = psf.from.env.entrez.fc(entrez.fc, kegg.collection, split)
  if(calculate.significance){
    cat("\nPerforming bootstrap calculations with", bst.steps, " steps\n")
    psf.results.sig = bootstrap.significance(psf.results, entrez.fc, bst.steps)
  }
  psf.resulsts.processed = process.psf.results(psf.results.sig)
  cat("PSF calculations finished!\n")
  return(psf.resulsts.processed)
}
