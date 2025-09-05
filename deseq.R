library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

setwd('/Users/ishitajain/Desktop/')
directory = '/Users/ishitajain/Documents/PhD_work/data_rna_2025/feature_counts/counts_file'
targets <- read.delim('combined_count_matrix.tsv',  header = TRUE, sep = '\t')

##============================================MAKING THE METADATA TABLE ====================================================================================

sample_files <- grep("\\.tsv$", list.files(directory), value = TRUE)
targets$group = factor(targets$group)
targets$donor = factor(targets$donor)
sampleTable <- data.frame(
  sampleName = targets$sample,
  fileName   = targets$file,
  condition  = targets$group,
  donor      = targets$donor,
  run        = targets$Run
)
rownames(sampleTable) <- sampleTable$sampleName


ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)

##================================================= MKAING THE BIO GROUP WHICH WILL BE USED FOR COLLAPSING THE TECHNICAL REPLICTAES ===========================================

sampleTable$bio_group <- paste(sampleTable$donor, sampleTable$condition, sep="_")


ddsHTSeq$bio_group <- factor(sampleTable$bio_group)


ddsHTSeq$condition <- factor(sampleTable$condition)

ddscoll = collapseReplicates(
  ddsHTSeq,
  groupby = ddsHTSeq$bio_group, 
  run = ddsHTSeq$run           
)

##===========================================================================================================================================================================
##Checking the collapsing works or not
cat('Sanity check if the collapsing works fine or not (OPTIONAL)')

collapsed_samples <- data.frame(sample = colnames(ddscoll),
                                bio_group = ddscoll$bio_group,
                                condition = ddscoll$condition,
                                stringsAsFactors = FALSE)
##====================================================================================================================================================================
cat('Samples before collapsing i.e. in ddsHTSeq \n')
print(ncol(ddsHTSeq))
cat('\n Samples before collapsing i.e. in ddscoll \n')
print(ncol(ddscoll))
cat('The number of samples belonging to the condition before collapsing \n')
print(table(ddsHTSeq$bio_group))
cat('The number of samples belonging to the condition before collapsing \n')
print(table(ddscoll$bio_group))
##=================================================================================== TO CHECK WHICH HOW THE TECHNICAL REPS HAVE COLLAPSED ================================================================================

original_samples <- data.frame(samples = colnames(ddsHTSeq),
                               condition = ddsHTSeq$condition,
                               bio_group = ddsHTSeq$bio_group,
                               run = ddsHTSeq$run,
                               stringsAsFactors = FALSE)



collapsed_samples <- data.frame(samples = colnames(ddscoll),
                               condition = ddscoll$condition,
                               bio_group = ddscoll$bio_group,
                               stringsAsFactors = FALSE)


for(group in unique(original_samples$bio_group)){
  sample_group <- original_samples[original_samples$bio_group == group,]
  if (nrow(sample_group) > 1){
    cat("\nBio_group '", group, "' collapsed from ", nrow(sample_group), " samples:\n", sep="")
    for (i in 1:nrow(sample_group)){
      cat("  - ", sample_group$sample[i], " (run: ",sample_group$run[i], ")\n", sep="")
    }
    cat("  → Into collapsed sample: ", group, "\n")
  }
}

##================================================================ CHECKING FOR CONTROL SAMPLES COLLAPSE SPECIFICALLY ================================================================

control_samples = original_samples[grepl('control', original_samples$bio_group, ignore.case = TRUE), ]
control_by_donor = table(control_samples$bio_group)

cat("\n=== VALIDATION SUMMARY ===\n")
cat("Total original samples:", ncol(ddsHTSeq), "\n")
cat("Final collapsed samples:", ncol(ddscoll), "\n")
cat("Expected final samples:", length(unique(original_samples$bio_group)), "\n")

if (ncol(ddscoll) == length(unique(original_samples$bio_group))){
  cat('=====Collapsing is Successful=====')
}else{
  cat('Please try again for collapsing as the number of conditions post and pre collapse do not match')
}

##=========================================================================================RUNNING DESEQ WITH COLLAPSED DATA ====================================================================================
ddscoll$condition <- relevel(ddscoll$condition, ref = "DSC_control")
dds <- DESeq(ddscoll)
res <- results(dds)
vsd <- vst(dds)
datExpr <- assay(vsd)
print(resultsNames(dds))
cond <- c("condition", "DSC_3x_UVB300_24h", "DSC_control")
  if(all(c('DSC_control',  'MZVH_control') %in% cond)){
  alpha = 0.01
}else{
  alpha = 0.05
}
  
res <- results(dds, contrast = cond)


upregulated_genes <- subset(res, res$padj<=alpha&res$log2FoldChange>0) 
downregulated_genes <- subset(res, res$padj<=alpha&res$log2FoldChange<0)
##===========================================================================================================================
print(length(rownames(upregulated_genes)))
print(length(rownames(downregulated_genes)))
####================================================ WRITING UP AND DPWN TO SEPERATE FILES ================================================
write.table(downregulated_genes, file='/Users/ishitajain/counts/downregulated_genes_DSC_vs_MZVH.tsv', sep = '\t')
write.table(upregulated_genes, file='/Users/ishitajain/counts/upregulated_genes_DSC_vs_MZVH.tsv', sep='\t')



