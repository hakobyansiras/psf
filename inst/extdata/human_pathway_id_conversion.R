kegg_human_pathway_list <- read.delim(file = "~/kegg_human_pathway_list.tsv", sep = "\t", stringsAsFactors = F)

kegg_human_pathway_list$Name <- gsub("(", "", kegg_human_pathway_list$Name, fixed = T)
kegg_human_pathway_list$Name <- gsub(")", "", kegg_human_pathway_list$Name, fixed = T)
kegg_human_pathway_list$Name <- gsub("-", "_", kegg_human_pathway_list$Name)
kegg_human_pathway_list$Name <- gsub(" ", "_", kegg_human_pathway_list$Name)
kegg_human_pathway_list$Name <- gsub("___", "_", kegg_human_pathway_list$Name)

pathway_name_to_id <- setNames(object = kegg_human_pathway_list$Entry, nm = kegg_human_pathway_list$Name)


save(pathway_name_to_id, file = "~/pathway_name_to_id.RData")
