#### Please download data file from Zenodo (https://zenodo.org/record/7396399/files/hsaSPIA.RData?download=1) and put it in R working directory
load(system.file("extdata", "melanoma_exp.RData", package="psf"))
load(system.file("extdata", "edited_pathways_new.RData", package="psf"))

### loading required libraries
library(DESeq2)
library(psf)
library(SPIA)
library(venn)

### running differential gene expression analysis with DESeq2
dds <- DESeqDataSetFromMatrix(round(melanoma_deseq_normalized_counts[,patient.concise$ID.indata[which(patient.concise$Group %in% c("Mel_1", "Naevi_1"))]]),
                              colData = patient.concise[which(patient.concise$Group %in% c("Mel_1", "Naevi_1")),],
                              design = ~ Group)
dds$Group <- factor(dds$Group, levels = c("Naevi_1","Mel_1"))
### setting size factors to 1 as the dataset is already DESeq normalized.
sizeFactors(dds) <- setNames(object = rep(1, 39), nm = patient.concise[which(patient.concise$Group %in% c("Mel_1", "Naevi_1")),"ID.indata"])
dds <- DESeq2::estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res_deseq <- results(dds)

## filtering significant genes
res_deseq_no_na <- res_deseq[which(!is.na(res_deseq$log2FoldChange)),]
significant_gene_fc <- setNames(nm = rownames(res_deseq_no_na)[which(res_deseq_no_na$padj < 0.05)], object = res_deseq_no_na$log2FoldChange[which(res_deseq_no_na$padj < 0.05)])

## filtering significant genes with high FC values
high_fc_genes <- res_deseq_no_na[which(abs(res_deseq_no_na$log2FoldChange) > 1),]
significant_gene_fc_high_fc <- setNames(nm = rownames(high_fc_genes)[which(high_fc_genes$padj < 0.05)], object = high_fc_genes$log2FoldChange[which(high_fc_genes$padj < 0.05)])

### SPIA with all significant genes. Please download hsaSPIA.RData file from zenodo and put it in working directory before ruunning this command.
res <- spia(de=significant_gene_fc, all = rownames(res_deseq_no_na), organism="hsa",data.dir="./")

res$pG=combfunc(res$pNDE,res$pPERT,combine="norminv")
res$pGFdr=p.adjust(res$pG,"fdr")
res$pGFWER=p.adjust(res$pG,"bonferroni")
plotP(res,threshold=0.05)

### SPIA with big FC significant genes. Please download hsaSPIA.RData file from zenodo and put it in working directory before ruunning this command.
res_high_fc <- spia(de=significant_gene_fc_high_fc, all = rownames(res_deseq_no_na), organism="hsa",data.dir="./")

res_high_fc$pG=combfunc(res_high_fc$pNDE, res_high_fc$pPERT,combine="norminv")
res_high_fc$pGFdr=p.adjust(res_high_fc$pG,"fdr")
res_high_fc$pGFWER=p.adjust(res_high_fc$pG,"bonferroni")
plotP(res_high_fc,threshold=0.05)

## runing PSF analysis with 2000 bst steps on the same pathway set with FC data
mel_psf <- run_psf(entrez.fc = matrix(data = 2^res_deseq_no_na$log2FoldChange, dimnames = list(rownames(res_deseq_no_na), "gene_fc")), 
                   kegg.collection = edited_pathways_new, calculate.significance = T, bst.steps = 2000, ncores = 4, shuffling_mode = "local")


## extracting significant pathway list
psf_sig_pathways <- unlist(lapply(mel_psf, function(x) {
  if(any(x$p_val_mat[x$sink.nodes,"gene_fc"] < 0.05)) {
    gsub("path:hsa", "", x$attr$name)
  }
}))


intersect(psf_sig_pathways, 
          res$ID[which(res$pGFdr < 0.05)])

l <- list(SPIA = res$ID[which(res$pGFdr < 0.05)], PSF = psf_sig_pathways)

l_high_fc_genes <- list(SPIA = res_high_fc$ID[which(res_high_fc$pGFdr < 0.05)], PSF = psf_sig_pathways)

### plotting venn diagramm
venn::venn(l, ggplot = F, zcolor = "style", box = F, ilabels = T, lty = 1, ilcs = 1.3, sncs = 1.5, lwd = 1.3)
venn::venn(l_high_fc_genes, ggplot = F, zcolor = "style", box = F, ilabels = T, lty = 1, ilcs = 1.3, sncs = 1.5, lwd = 1.3)


spia_specific_pathways <- Reduce(res$ID[which(res$pGFdr < 0.05)],
                                 psf_sig_pathways)

psf_specific_pathways <- setdiff(psf_sig_pathways, 
                                 res$ID[which(res$pGFdr < 0.05)])

### plotting AMPK signaling pathway
plot_pathway(mel_psf$AMPK_signaling_pathway, plot_type = "kegg", color_nodes = "psf_activities", use_old_images = T, log_norm = T)