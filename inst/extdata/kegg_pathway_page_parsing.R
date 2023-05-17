library(rvest)

# Load the webpage
webpage <- read_html("https://www.kegg.jp/kegg/pathway.html")

# Extract the headings and subheadings
headings_and_subheadings <- html_nodes(webpage, 'b') %>% html_text()

headings_and_subheadings <- headings_and_subheadings[grep("\\d+", headings_and_subheadings)]

# Extracting subheading content
subheadings_content <- webpage %>% html_nodes(".list") %>% html_text2()

subheadings_content <- lapply(subheadings_content, function(x) {
  
  text_vec <- unlist(strsplit(x, split = "\n"))
  
  t(sapply(split(text_vec, rep(1:ceiling(length(text_vec)/2), each=2, length.out=length(text_vec))), function(y) {
    if(length(unlist(strsplit(y[1], split = " "))) > 1) {
      c(unlist(strsplit(y[1], split = " "))[1],
        paste(unlist(strsplit(y[1], split = " "))[2:length(unlist(strsplit(y[1], split = " ")))], collapse = " "),
        y[2])
    } else {
      c(unlist(strsplit(y[1], split = " "))[1],
        "",
        y[2])
    }
    
  }))
  
})

### hardcoded here
headings <- headings_and_subheadings[1:7]
subheadings <- headings_and_subheadings[8:length(headings_and_subheadings)]

kegg_global_pathways_table <- Reduce(rbind, lapply(1:length(subheadings_content), function(x) {
  
  cbind(subheadings_content[[x]], 
        rep(subheadings[x], nrow(subheadings_content[[x]])), 
        rep(headings[grep(substr(subheadings[x], 1, 2), headings)], nrow(subheadings_content[[x]]) ))
  
}))

colnames(kegg_global_pathways_table) <- c("id", "avail_database", "name", "subclass", "class")
rownames(kegg_global_pathways_table) <- NULL
kegg_global_pathways_table <- kegg_global_pathways_table[,c("class", "subclass", "name", "id", "avail_database")]

write.table(kegg_global_pathways_table, file = "~/kegg_global_pathways_table.tsv", sep = "\t", quote = F, row.names = F)
