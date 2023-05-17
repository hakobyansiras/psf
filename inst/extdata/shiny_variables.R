kegg_human_pathway_list <- read.delim(file = "~/kegg_human_pathway_list.tsv", sep = "\t", stringsAsFactors = F)

kegg_human_pathway_list$Name <- gsub("(", "", kegg_human_pathway_list$Name, fixed = T)
kegg_human_pathway_list$Name <- gsub(")", "", kegg_human_pathway_list$Name, fixed = T)
kegg_human_pathway_list$Name <- gsub("-", "_", kegg_human_pathway_list$Name)
kegg_human_pathway_list$Name <- gsub(" ", "_", kegg_human_pathway_list$Name)
kegg_human_pathway_list$Name <- gsub("___", "_", kegg_human_pathway_list$Name)
kegg_human_pathway_list$Name <- gsub("/", "", kegg_human_pathway_list$Name)
kegg_human_pathway_list$Name <- gsub("__", "_", kegg_human_pathway_list$Name)

pathway_codes_new <- setNames(object = kegg_human_pathway_list$Entry, nm = kegg_human_pathway_list$Name)

kegg_arrows_type <- c("arrow", "arrow", "arrow", "arrow", "bar", "arrow", "arrow", "arrow", "bar", "arrow", "", "arrow", "arrow", "bar", "arrow", "arrow")
names(kegg_arrows_type) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
line_col <- c("red", "red", "red", "red", "blue","red", "red", "red", "blue", "red", "red", "red", "red", "blue", "red", "red")
names(line_col) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )


save(all_interactions, color_legend_maker, edge_subtype1, edge_subtype2, entrez_to_symbol, kegg_compounds_to_full_name, pal1, pal2, 
     pathway_codes_new, kegg_human_pathway_list, line_col, kegg_arrows_type, file = "inst/shinyApp/whole_data_unit.RData")

