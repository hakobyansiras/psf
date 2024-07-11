library(data.table)

signor_all_pathways <- fread("~/Signor_database/all_signor_pathways.tsv", sep = "\t", stringsAsFactors = F, data.table = F)

signor_all_pathways$pathway_name <- gsub(" ", "_", signor_all_pathways$pathway_name)

signor_pathwa_names <- unique(signor_all_pathways$pathway_name)

# library(biomaRt)
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# gene_symbol_to_entrez <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), mart = ensembl)
# gene_symbol_to_entrez <- gene_symbol_to_entrez[which(!is.na(gene_symbol_to_entrez$entrezgene_id)),]
# gene_symbol_to_entrez <- gene_symbol_to_entrez[which(!duplicated(gene_symbol_to_entrez$entrezgene_id)),]
# gene_symbol_to_entrez <- setNames(object = gene_symbol_to_entrez$entrezgene_id, nm = gene_symbol_to_entrez$hgnc_symbol)
# 
# uniprot_to_entrez <- getBM(attributes=c('entrezgene_id', 'uniprot_gn_id'), mart = ensembl)
# uniprot_to_entrez <- uniprot_to_entrez[which(!is.na(uniprot_to_entrez$entrezgene_id)),]
# uniprot_to_entrez <- uniprot_to_entrez[which(!duplicated(uniprot_to_entrez$entrezgene_id)),]
# uniprot_to_entrez <- setNames(object = uniprot_to_entrez$entrezgene_id, nm = uniprot_to_entrez$uniprot_gn_id)

signor_type_converter <- setNames(nm = c("antibody", "chemical", "complex", "fusion protein", "mirna", "phenotype", "protein", "proteinfamily", "smallmolecule", "stimulus"),
                                  object = c("gene", "compound", "gene", "gene", "compound", "phenotype", "gene", "gene", "compound", "stimulus"))


signor_PFs <- read.delim(file = "~/Signor_database/SIGNOR_PF.csv", sep = ";", stringsAsFactors = F)
signor_PFs_vec <- setNames(object = signor_PFs$LIST.OF.ENTITIES, nm = signor_PFs$SIGNOR.ID)

signor_complexes <- read.delim(file = "~/Signor_database/SIGNOR_complexes.csv", sep = ";", stringsAsFactors = F)
signor_complexes_vec <- setNames(object = signor_complexes$LIST.OF.ENTITIES, nm = signor_complexes$SIGNOR.ID)


signor_all_pathways_processed <- signor_all_pathways

signor_all_pathways_processed$ida[which(signor_all_pathways_processed$typea == "proteinfamily")] <-  unname(signor_PFs_vec[signor_all_pathways_processed$ida[which(signor_all_pathways_processed$typea == "proteinfamily")]])
signor_all_pathways_processed$idb[which(signor_all_pathways_processed$typeb == "proteinfamily")] <-  unname(signor_PFs_vec[signor_all_pathways_processed$idb[which(signor_all_pathways_processed$typeb == "proteinfamily")]])

signor_all_pathways_processed$ida[which(signor_all_pathways_processed$typea == "complex")] <-  unname(signor_complexes_vec[signor_all_pathways_processed$ida[which(signor_all_pathways_processed$typea == "complex")]])
signor_all_pathways_processed$idb[which(signor_all_pathways_processed$typeb == "complex")] <-  unname(signor_complexes_vec[signor_all_pathways_processed$idb[which(signor_all_pathways_processed$typeb == "complex")]])


signor_uniprot_ids <- unique(
  c(
    unlist(lapply(signor_all_pathways_processed$ida[which(signor_all_pathways_processed$typea %in% c("protein", "proteinfamily", "complex"))], function(x) {
      unlist(strsplit(x, split = ",  "))
    })),
    unlist(lapply(signor_all_pathways_processed$idb[which(signor_all_pathways_processed$typeb %in% c("protein", "proteinfamily", "complex"))], function(x) {
      unlist(strsplit(x, split = ",  "))
    }))
  )
)

## writing signor uniprot ids to use web based mapping tool for uniprot to entrez conversion
write.table(signor_uniprot_ids, file = "~/signor_uniprot_ids.txt", quote = F, sep = "", row.names = F)
uniprot_to_entrez_new <- read.delim(file = "~/idmapping_2023_09_22.tsv", sep = "\t", stringsAsFactors = F)
uniprot_to_entrez_new <- setNames(object = uniprot_to_entrez_new$To, nm = uniprot_to_entrez_new$From)


signor_all_pathways_processed$ida[which(signor_all_pathways_processed$typea %in% c("protein", "proteinfamily", "complex"))] <- unname(sapply(signor_all_pathways_processed$ida[which(signor_all_pathways_processed$typea %in% c("protein", "proteinfamily", "complex"))], function(x) {
  gene_ids <- uniprot_to_entrez_new[unlist(strsplit(x, split = ",  "))]
  gene_ids[which(is.na(gene_ids))] <- 0
  paste(gene_ids, collapse = ",")
}))

signor_all_pathways_processed$idb[which(signor_all_pathways_processed$typeb %in% c("protein", "proteinfamily", "complex"))] <- unname(sapply(signor_all_pathways_processed$idb[which(signor_all_pathways_processed$typeb %in% c("protein", "proteinfamily", "complex"))], function(x) {
  gene_ids <- uniprot_to_entrez_new[unlist(strsplit(x, split = ",  "))]
  gene_ids[which(is.na(gene_ids))] <- 0
  paste(gene_ids, collapse = ",")
}))

interaction_type_coverter <- setNames(nm = c("up-regulates quantity", "down-regulates quantity", "up-regulates", "up-regulates quantity by expression", "down-regulates quantity by repression", "down-regulates", "up-regulates activity", "down-regulates quantity by destabilization", "down-regulates activity", "unknown", "form complex", "up-regulates quantity by stabilization"), 
                                      object = c("activation", "inhibition", "activation", "activation", "inhibition", "inhibition", "activation", "inhibition", "inhibition", "activation", "activation", "activation"))

signor_all_pathways_processed$type <- interaction_type_coverter[signor_all_pathways_processed$effect]

signor_all_pathways_processed$type_a <- signor_type_converter[signor_all_pathways_processed$typea]
signor_all_pathways_processed$type_b <- signor_type_converter[signor_all_pathways_processed$typeb]


signor_pathway_import <- function(pathway_table) {
  
  pathway_table <- pathway_table[which(pathway_table$entitya != pathway_table$entityb),]
  pathway_table <- pathway_table[which(!duplicated(pathway_table[,c("entitya", "entityb")])),]
  
  table_a <- pathway_table[,c("entitya", "type_a", "ida", "databasea")] 
  colnames(table_a) <- c("label", "type", "signor_id", "database")
  
  table_b <- pathway_table[,c("entityb", "type_b", "idb", "databaseb")] 
  colnames(table_b) <- c("label", "type", "signor_id", "database")
  
  node_components <- rbind(table_a, table_b)
  
  node_components <- node_components[which(!duplicated(node_components[,1:2])),]

  colnames(node_components) <- c("label", "type", "component_id_s", "database")
  
  node_table <- data.frame(node_id = as.character(1:nrow(node_components)),
                           label = node_components$label,
                           component_id_s = node_components$component_id_s,
                           type = node_components$type,
                           psf_function = "mean",
                           expression = 1,
                           signal = 1,
                           # x = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.x"]))),
                           # y = as.integer(unlist(graph::nodeData(pathway$graph, attr = attribute_converter["kegg.gr.y"]))),
                           existence = "exist",
                           change_info = "no_change",
                           data_source = "signor",
                           stringsAsFactors = F
  )
  
  label_to_id <- setNames(nm = node_table$label, object = node_table$node_id)
  
  edge_table <- data.frame(from = label_to_id[pathway_table$entitya],
                           to = label_to_id[pathway_table$entityb],
                           type = "",
                           subtype1 = pathway_table$type,
                           subtype2 = ifelse(pathway_table$direct == "t", "indirect effect", "direct"),
                           # subtype2 = unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype2")),
                           state = "",
                           weight = 1,
                           stringsAsFactors = F
  )
  
  return(list(node_table = node_table, edge_table = edge_table))
  
}


signor_collection <- lapply(signor_pathwa_names, function(x) {
  
  print(x)
  
  signor_tables <- signor_pathway_import(pathway_table = signor_all_pathways_processed[which(signor_all_pathways_processed$pathway_name == x),])
  
  pathway <- df_to_graphnel(node_table = signor_tables$node_table, edge_table = signor_tables$edge_table)
  
  pathway$attrs$title <- x
  
  pathway
  
})

names(signor_collection) <- signor_pathwa_names

### removing Glycogenesis pathway as it does not have sink nodes
signor_collection <- signor_collection[names(signor_collection)[which(names(signor_collection) != "Glycogenesis")]]

