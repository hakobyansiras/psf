#### Please download data files from Zenodo () and put the files in R working directory
#### NOTE you can load preprocessed big files from RData file with command bellow
load("isoform_analysis_files.RData")

### loading required libraries
library(tximport)
library(DESeq2)
library(psf)
library(biomaRt)
library(rtracklayer)
library(IsoformSwitchAnalyzeR)
library(venn)
library(plotly)
library(ggplot2)

### loading annotation (unccomment and run if you prefer to run preprocessing too)
# gencode32 <- rtracklayer::import('../isoform_bmu_calculation_for_melanoma/gencode.v32.chr_patch_hapl_scaff.annotation.gtf')
# 
# gencode32 <- rtracklayer::as.data.frame(gencode32)
# 
# gencode32 <- gencode32[which(gencode32$type == "transcript"),]
# 
# gencode32 <- gencode32[grep("chr", gencode32$seqnames),]

transcript_to_ensembl_gene <- setNames(object = gencode32$gene_id, nm = gencode32$transcript_id)

### importing kallisto abundance qunatification data with tximport library (unccomment and run if you prefer to run preprocessing too)
# files = file.path("melanoma_new_quantification", list.files("melanoma_new_quantification"), "quantification", "abundance.tsv")
# names(files) = list.files('melanoma_new_quantification')
# 
# tx2gene <- data.frame(isoform = gencode32$transcript_id, gene = gencode32$gene_id)
# 
# isoform_aggregated_list <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

## Sample name corrections (unccomment and run if you prefer to run preprocessing too)
# samples = gsub("_", ".", colnames(isoform_aggregated_list$counts)) 
# samples[ grepl("^[0-9]", samples)] = paste0("X", samples[ grepl("^[0-9]", samples)])
# 
# colnames(isoform_aggregated_list$counts) <- samples
# colnames(isoform_aggregated_list$abundance) <- samples
# colnames(isoform_aggregated_list$length) <- samples


### loading patient data
load("patient.concise.RData")

### Normalizing expression count data with DESeq2 package
## deseq normalization of aggregated gene counts
dds <- DESeqDataSetFromTximport(isoform_aggregated_list,
                                colData = patient.concise,
                                design = ~ 1)

dds <- estimateSizeFactors(dds)

normalized_aggregated_counts <- counts(dds, normalized=TRUE)

rownames(normalized_aggregated_counts) <- substr(rownames(normalized_aggregated_counts), 1, 15)

### Calculating coding potatntial of isofroms from transcripts sequence files (unccomment and run if you prefer to run preprocessing too)
# system("/storage/users/siras/Isoform_psf_calculation/CPC2/bin/CPC2.py -i gencode.v32.transcripts.fa -o coding_gencod32_transcript.txt")
# 
# cpc2_coding_transcripts <- read.delim(file = "coding_gencod32_transcript.txt", sep = "\t", stringsAsFactors = F)
# 
# cpc2_coding_transcripts$transcript_id <- sapply(cpc2_coding_transcripts$X.ID, function(x) {unlist(strsplit(x, split = "|", fixed = T))[1]})


### Importing isoform expression data (unccomment and run if you prefer to run preprocessing too)
# kallistoQuant <- importIsoformExpression(
#   sampleVector = files,
#   addIsofomIdAsColumn = TRUE
# )

### Calculating isoform fractions
isoform_data_IF <- as.data.frame(kallistoQuant$abundance)

isoform_data_IF <- cbind(isoform_data_IF$isoform_id, transcript_to_ensembl_gene[isoform_data_IF$isoform_id], isoform_data_IF[,-1])

colnames(isoform_data_IF)[1:2] <- c("isoform_id", "gene_id")

rownames(isoform_data_IF) <- NULL

isoform_data_IF$isoform_id <- as.character(isoform_data_IF$isoform_id)

isoform_data_IF$gene_id <- as.character(isoform_data_IF$gene_id)

isoform_fractions <- isoformToIsoformFraction(isoformRepExpression = isoform_data_IF)

rownames(isoform_fractions) <- isoform_fractions$isoform_id

isoform_fractions <- as.matrix(isoform_fractions[,-1])

isoform_fractions[which(is.nan(isoform_fractions), arr.ind = T)] <- 0

### fraction extraction and exp procesing function

isoform_fraction_extractor <- function(isoform_fractions, aggregated_counts, non_coding_transcripts, ensembl_to_entrez, transcript_to_ensembl_gene) {
  
  transcript_to_ensembl_gene <- substr(transcript_to_ensembl_gene, 1, 15)
  
  non_coding_isoform_fractions <- isoform_fractions[non_coding_transcripts,]
  
  non_coding_isoform_fractions_aggregated <- aggregate(x = non_coding_isoform_fractions, by = list(transcript_to_ensembl_gene[rownames(non_coding_isoform_fractions)]), FUN = sum)
  
  rownames(non_coding_isoform_fractions_aggregated) <- non_coding_isoform_fractions_aggregated$Group.1
  
  non_coding_isoform_fractions_aggregated <- as.matrix(non_coding_isoform_fractions_aggregated[,-1])
  
  non_coding_isoform_fractions_aggregated <- 1 - non_coding_isoform_fractions_aggregated
  
  non_coding_isoform_fractions_aggregated[which(non_coding_isoform_fractions_aggregated < 0)] <- 0
  
  non_coding_isoform_fractions_aggregated <- non_coding_isoform_fractions_aggregated + 0.00001
  
  non_coding_isoform_fractions_aggregated[which(non_coding_isoform_fractions_aggregated > 1)] <- 1
  
  ### fraction extraction
  non_coding_isoform_fraction_extracted_gene_exp <- aggregated_counts
  
  non_coding_isoform_fraction_extracted_gene_exp[rownames(non_coding_isoform_fractions_aggregated),] <- non_coding_isoform_fraction_extracted_gene_exp[rownames(non_coding_isoform_fractions_aggregated),]*non_coding_isoform_fractions_aggregated
  
  
  ## keeping genes which have been expressed at least in 2 samples
  normalized_aggregated_counts_filtered <- aggregated_counts[unname(unlist(sapply(rownames(aggregated_counts), function(x) {if(sum(aggregated_counts[x,] > 0) >= 2) {x}}))),]
  non_coding_isoform_fraction_extracted_gene_exp <- non_coding_isoform_fraction_extracted_gene_exp[unname(unlist(sapply(rownames(non_coding_isoform_fraction_extracted_gene_exp), function(x) {if(sum(non_coding_isoform_fraction_extracted_gene_exp[x,] > 0) >= 2) {x}}))), ]
  
  
  ensembl_to_entrez_filtered <- ensembl_to_entrez[which(ensembl_to_entrez$ensembl_gene_id %in% rownames(normalized_aggregated_counts_filtered)),]
  
  ensembl_to_entrez_filtered <- ensembl_to_entrez_filtered[which(!duplicated(ensembl_to_entrez_filtered$entrezgene_id)),]
  
  ensembl_to_entrez_filtered <- setNames(object = ensembl_to_entrez_filtered$entrezgene_id, nm = ensembl_to_entrez_filtered$ensembl_gene_id)
  
  
  normalized_aggregated_counts_filtered <- normalized_aggregated_counts_filtered[names(ensembl_to_entrez_filtered),]
  non_coding_isoform_fraction_extracted_gene_exp <- non_coding_isoform_fraction_extracted_gene_exp[names(ensembl_to_entrez_filtered),]
  
  
  rownames(normalized_aggregated_counts_filtered) <- ensembl_to_entrez_filtered[rownames(normalized_aggregated_counts_filtered)]
  rownames(non_coding_isoform_fraction_extracted_gene_exp) <- ensembl_to_entrez_filtered[rownames(non_coding_isoform_fraction_extracted_gene_exp)]
  
  return(list(gene_exp = normalized_aggregated_counts_filtered, isoform_exp = non_coding_isoform_fraction_extracted_gene_exp,
              aggregated_frations = non_coding_isoform_fractions_aggregated))
}

### Downloading ensembl to entrez conversion list
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensembl_to_entrez <- getBM(attributes=c('entrezgene_id', 'ensembl_gene_id'), 
                               mart = ensembl)
ensembl_to_entrez <- ensembl_to_entrez[which(!is.na(ensembl_to_entrez$entrezgene_id)), ]

### calculating coding gene expression fractions
cpc_fractions <- isoform_fraction_extractor(isoform_fractions = isoform_fractions, aggregated_counts = normalized_aggregated_counts, 
                                            non_coding_transcripts = cpc2_coding_transcripts$transcript_id[which(cpc2_coding_transcripts$label == "noncoding")], 
                                            ensembl_to_entrez = ensembl_to_entrez, transcript_to_ensembl_gene = transcript_to_ensembl_gene)


##### stopped here remove NA values from expression data, may be include it into above function

### Differential gene expression analysis between Melanoma type 1 and Naevi type 1
dds_gene <- DESeqDataSetFromMatrix(round(cpc_fractions$gene_exp[,patient.concise$ID.indata[which(patient.concise$Group %in% c("Mel_1", "Naevi_1"))]]),
                                   colData = patient.concise[which(patient.concise$Group %in% c("Mel_1", "Naevi_1")),],
                                   design = ~ Group)
dds_gene$Group <- factor(dds_gene$Group, levels = c("Naevi_1","Mel_1"))
sizeFactors(dds_gene) <- setNames(object = rep(1, 39), nm = patient.concise[which(patient.concise$Group %in% c("Mel_1", "Naevi_1")),"ID.indata"])
dds_gene <- DESeq2::estimateDispersions(dds_gene)
dds_gene <- nbinomWaldTest(dds_gene)
res_gene <- results(dds_gene)


dds_isoform <- DESeqDataSetFromMatrix(round(cpc_fractions$isoform_exp[,patient.concise$ID.indata[which(patient.concise$Group %in% c("Mel_1", "Naevi_1"))]]),
                                      colData = patient.concise[which(patient.concise$Group %in% c("Mel_1", "Naevi_1")),],
                                      design = ~ Group)
dds_isoform$Group <- factor(dds_isoform$Group, levels = c("Naevi_1","Mel_1"))
sizeFactors(dds_isoform) <- setNames(object = rep(1, 39), nm = patient.concise[which(patient.concise$Group %in% c("Mel_1", "Naevi_1")),"ID.indata"])
dds_isoform <- DESeq2::estimateDispersions(dds_isoform)
dds_isoform <- nbinomWaldTest(dds_isoform)
res_isoform <- results(dds_isoform)

venn::venn(list("Full gene" = rownames(res_gene)[which(res_gene$padj < 0.05)],
                "Coding gene" = rownames(res_isoform)[which(res_isoform$padj < 0.05)]), 
           ggplot = F, zcolor = "style", box = F, ilabels = T, lty = 1, ilcs = 1.3, sncs = 1.5, lwd = 1.3)


isoform_fc_matrix <- matrix(data = 2^res_isoform$log2FoldChange, dimnames = list(rownames(res_isoform), "isoform_fc"))
isoform_fc_matrix <- isoform_fc_matrix[which(!is.na(isoform_fc_matrix[,1])),, drop = F]

gene_fc_matrix <- matrix(data = 2^res_gene$log2FoldChange, dimnames = list(rownames(res_gene), "gene_fc"))
gene_fc_matrix <- gene_fc_matrix[which(!is.na(gene_fc_matrix[,1])),, drop = F]


## loading curated pathway collection
load(system.file("extdata", "edited_pathways_new.RData", package="psf"))


#### calculating psf for full and coding gene expression Deseq FC data (calculation will take rougly 8 hour, you can load precalculate data with command below) ####
load("psf_with_significance.RData")

## uncomment the script to run psf analysis
# load("psf_significance_calc_functions.RData")
#
# mel_psf_single_norm <- psf_calc_with_significance(fc_matrix = gene_fc_matrix, 
#                                                   pathway_collection = edited_pathways_new, 
#                                                   shuffle_mode = "all", bst_steps = 2000)
# 
# mel_fraction_extracted_psf_single_norm <- psf_calc_with_significance(fc_matrix = isoform_fc_matrix, 
#                                                                      pathway_collection = edited_pathways_new, 
#                                                                      shuffle_mode = "all", bst_steps = 2000)



#### Fulll and coding gene expression psf results comparison ####
sink_labels <- Reduce(c,
                      lapply(names(mel_psf_single_norm$psf_list$gene_fc), function(x) {
                        pathway_mat <- sapply(mel_psf_single_norm$psf_list, function(y) {
                          y[[x]]$signal.at.sink
                        })
                        paste(x, rownames(pathway_mat), sep = "_")
                      }))


psf_data_gene <- mel_psf_single_norm$psf_matrix

psf_data_gene <- psf_data_gene[which(psf_data_gene$node_label %in% sink_labels),]

rownames(psf_data_gene) <- psf_data_gene$node_label

psf_data_gene_filtered <- psf_data_gene[which(psf_data_gene$gene_fc_pval < 0.05), ]

psf_data_isoform <- mel_fraction_extracted_psf_single_norm$psf_matrix

psf_data_isoform <- psf_data_isoform[which(psf_data_isoform$node_label %in% sink_labels),]

rownames(psf_data_isoform) <- psf_data_isoform$node_label

psf_data_isoform_filtered <- psf_data_isoform[which(psf_data_isoform$isoform_fc_pval < 0.05), ]

### building venn diagram to comaper full and coding gene expression based significant pathway sink nodes
venn::venn(list("Full gene" = psf_data_gene_filtered$node_label, "Coding gene" = psf_data_isoform_filtered$node_label), 
           ggplot = F, zcolor = "style", box = F, ilabels = T, lty = 1, ilcs = 1.3, sncs = 1.5, lwd = 1.3)

### Scatterplot for full and coding gene expression PSF compariosn
ggplot(data.frame(Full_gene = log10(psf_data_gene$gene_fc), Coding_gene = log10(psf_data_isoform$isoform_fc)), aes(x=Full_gene, y=Coding_gene)) + 
  geom_point(size = 3) + xlab("Full gene log10 PSF values") + ylab("Coding gene log10 PSF values") + 
  theme(axis.title = element_text(size = 32), axis.text = element_text(size = 28))


### filtering only coding expression significant pathway sinks
isoform_specific_dereg_sinks <- setdiff(psf_data_isoform_filtered$node_label, psf_data_gene_filtered$node_label)
gene_specific_dereg_sinks <- setdiff(psf_data_gene_filtered$node_label, psf_data_isoform_filtered$node_label)

full_vs_coding_psf <- cbind(psf_data_isoform_filtered[isoform_specific_dereg_sinks,], psf_data_gene[isoform_specific_dereg_sinks, 2:3],  psf_data_isoform_filtered[isoform_specific_dereg_sinks, "isoform_fc"] - psf_data_gene[isoform_specific_dereg_sinks, "gene_fc"])

colnames(full_vs_coding_psf) <- c("node_label", "coding_exp_psf", "coding_exp_psf_pval", "full_exp_psf", "full_exp_psf_pval", "psf_diff")

View(full_vs_coding_psf)

#### pathway visualization function
pathway_vis <- function(pathway_name, deseq_deg_list, psf_list, sig_sink_list, plot_type, log_norm) {
  
  sig_degs <- rownames(deseq_deg_list)[which(deseq_deg_list$padj < 0.05)]
  
  sig_deg_nodes <- unname(unlist(sapply(graph::nodes(psf_list$psf_list[[1]][[pathway_name]]$graph), function(x) {
    if(any(unlist(graph::nodeData(psf_list$psf_list[[1]][[pathway_name]]$graph, x, attr = "genes")) %in% sig_degs)) {
      x
    }
  })))
  
  significant_sink_ids <- sapply(sig_sink_list[grep(pathway_name, sig_sink_list)], function(x) {
    unlist(strsplit(x, split = "_"))[length(unlist(strsplit(x, split = "_")))]
  })
  
  bot_sig_nodes <- intersect(sig_deg_nodes, significant_sink_ids)
  
  only_degs <- setdiff(sig_deg_nodes, bot_sig_nodes)
  only_sinks <- setdiff(significant_sink_ids, bot_sig_nodes)
  
  highlight_nodes <- c(only_degs, only_sinks, bot_sig_nodes)
  
  higlight_colors <- c(rep("red", length(only_degs)), rep("green", length(only_sinks)), rep("yellow", length(bot_sig_nodes)))
  
  psf::plot_kegg_image_pathway(psf_list$psf_list[[1]][[pathway_name]], no_color_mode = F, 
                               mapping_data_type = "signal", use_old_images = T, 
                               highlight_nodes = highlight_nodes, highlight_color = higlight_colors,
                               y_adj_text = 13, y_adj_sink = 38, plot_type = plot_type, log_norm = log_norm
  )
  
}


## red - differentially expressed significant genes according to coding gene expression deseq
## green significantly deregulated PSF sink nodes
## yellow both PSF and DESeq significant
### plotting only coding gene expression significant sinks on cAMP_signaling_pathway pathway
pathway_vis(pathway_name = "cAMP_signaling_pathway", deseq_deg_list = res_isoform, 
            psf_list = mel_fraction_extracted_psf_single_norm, sig_sink_list = isoform_specific_dereg_sinks, plot_type = "kegg", log_norm = F)

### plotting only full gene expression significant sinks on FoxO_signaling_pathway signaling pathway
pathway_vis(pathway_name = "FoxO_signaling_pathway", deseq_deg_list = res_gene, 
            psf_list = mel_psf_single_norm, sig_sink_list = gene_specific_dereg_sinks, plot_type = "kegg", log_norm = F)



### PSF activity difference barplot builder
psf_diff_barplot_builder <- function(pathway_name) {
  
  gene_and_isoform_psf_mat <- cbind(psf_data_isoform_filtered[isoform_specific_dereg_sinks[grep(pathway_name, isoform_specific_dereg_sinks)],], 
                                    psf_data_gene[isoform_specific_dereg_sinks[grep(pathway_name, isoform_specific_dereg_sinks)],]
  )
  
  gene_and_isoform_psf_mat$isoform_fc <- log(gene_and_isoform_psf_mat$isoform_fc)
  gene_and_isoform_psf_mat$gene_fc <- log(gene_and_isoform_psf_mat$gene_fc)
  
  
  gene_and_isoform_psf_mat$node_label <- unlist(sapply(gene_and_isoform_psf_mat$node_label, function(x) {
    node_id <- unlist(strsplit(x, split = "_"))[length(unlist(strsplit(x, split = "_")))]
    
    graph::nodeData(edited_pathways_new[[pathway_name]]$graph, node_id, attr = "label")
    
  }))
  
  library(plotly)
  plot_ly(gene_and_isoform_psf_mat, x = ~node_label, y = ~gene_fc, type = 'bar', name = 'Full exp PSF', marker = list(color = "#0099ff")) %>%
    add_trace(y = ~isoform_fc, name = 'Coding exp PSF', marker = list(color = "#003399")) %>% 
    layout(yaxis = list(title = 'Log PSF value'), xaxis = list(title = 'Sink node'), barmode = 'group', title = pathway_name)
  
}

### plotting psf differences
psf_diff_barplot_builder("cAMP_signaling_pathway")
