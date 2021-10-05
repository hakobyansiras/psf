selected_pathway_set <- read.delim(file = "/home/siras/PSFC/inst/extdata/kgmls/selected_pathways.txt", sep = "\t", stringsAsFactors = F)


selected_pathway_set <- data.frame(pathway_code = selected_pathway_set$pathway_code,
                                   kgml_link = paste0("http://rest.kegg.jp/get/", selected_pathway_set$pathway_code, "/kgml"),
                                   png_link = paste0("http://rest.kegg.jp/get/", selected_pathway_set$pathway_code, "/image"),
                                   stringsAsFactors = F
                                   )


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


psf::generate.kegg.collection.from.kgml()