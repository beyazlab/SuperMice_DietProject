export LC_ALL=C
module load EBModules
module load R/4.1.0-foss-2021a

R
#devtools::install_github('kevinblighe/PCAtools')
#devtools::install_github("mkearney/rmd2jupyter")

library(circlize)

library('edgeR')
library("DESeq2")
#cpm = edgeR::cpm
library('PCAtools')
biplot = PCAtools::biplot
library(biomaRt)
library(tibble)
library(ggplot2)
library(gplots)
library(patchwork)
library(gridExtra)
library(ComplexHeatmap)
library(data.table)
library("grid")
library(ggrepel)
library("cowplot")
library("ggpubr")
library(dplyr)
library("apeglm")
library(UpSetR)
library("vsn")

proj_path <- "/grid/beyaz/home/subhash/projects/Diet_cohorts/BulkRNAseq/"
inpath <- paste0(proj_path,"quantification/featureCounts/")
outpath <- paste0(proj_path,"DE_results/")
plot_path <- paste0(outpath,"plots/")
frich_path <- paste0(outpath,"functional_analysis/")
meta_path <- "/grid/beyaz/home/subhash/projects/Diet_cohorts/metadata/"

x <- read.table(paste0(inpath,"merged_gene_counts_matchingRRBS_clean.txt"), header=T, sep="\t", row.names="Geneid")
xanot <- read.table(paste0(inpath,"merged_gene_counts_matchingRRBS_clean.txt"), header=T, sep="\t")
mdata <- read.table(paste0(meta_path,"supermice_dietcohort_complete_metadata_RRBS_RNAseq.txt"), header=T, sep="\t")


x <- x[,2:ncol(x)]


pdf(paste0(plot_path,"Diet_cohort_bulkRNAseq_sample_stats.pdf"), height=10); 
mat <- as.data.frame(mdata)
mat$tissue <- factor(mat$tissue)
p1 <- ggplot(mat, aes(x=current_Diet)) + geom_bar(stat="count", col="#2E2E2E") + geom_text(aes(label = ..count..), stat = "count", position = position_stack(vjust = 0.5), col="white") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
p2 <- ggplot(mat, aes(x=reversal_diet, fill = Sex, group=Sex)) + geom_bar(stat="count") + geom_text(aes(label = ..count..), stat = "count", position = position_stack(vjust = 0.5)) +
scale_fill_manual(name="Sex",values=c("#FFBF00","#8904B1")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~duration_of_reversion_months) 
p3 <- ggplot(mat, aes(x=reversal_diet, fill = Sex, group=Sex)) + geom_bar(stat="count") + geom_text(aes(label = ..count..), stat = "count", position = position_stack(vjust = 0.5)) +
scale_fill_manual(name="Sex",values=c("#FFBF00","#8904B1")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~duration_of_reversion_group) 

p4 <- ggplot(mat, aes(x=reversal_status, fill = Sex, group=Sex)) + geom_bar(stat="count") + geom_text(aes(label = ..count..), stat = "count", position = position_stack(vjust = 0.5)) +
scale_fill_manual(name="Sex",values=c("#FFBF00","#8904B1")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~original_Diet)


p5 <- ggplot(mat, aes(original_Diet, fill = factor(reversal_status)))+
  geom_bar(stat="count") +
  geom_text(aes(label = ..count..), stat = "count", position = "fill", vjust=3) +
  scale_fill_manual(name="Reversal",values=c("#FFBF00","#8904B1")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~duration_of_reversion_group)
  
  ggarrange(p2, p5, nrow=2, ncol=1)
  ggarrange(p1, p2, nrow=2, ncol=1)
  ggarrange(p3, p4, nrow=2, ncol=1)
  
dev.off()










####################### BOTTOM CODE NOT USED 

x <- x[,2:ncol(x)]
mdata <- as.data.frame(mdata)
rownames(mdata) <- mdata$unique_sample_ids
mdata <- mdata[, -which(names(mdata) == "sample_notes")]
mdatax <- subset(mdata,mdata$original_Diet %in% c("control","coconut"))
xm <- x[,mdatax$unique_sample_ids]
mdatax <- mdatax[match(colnames(xm), mdatax$unique_sample_ids),]
table(rownames(mdatax)==colnames(xm))
mdatax$group <- factor(paste0(mdatax$original_Diet,"_",mdatax$duration_of_reversion_group))

dds <- DESeqDataSetFromMatrix(countData = xm, colData = mdatax, design = ~Sex+group)

featureData <- data.frame(gene=rownames(xm))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)			
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


dds$condition <- factor(dds$group, levels = levels(dds$group))
#dds$condition <- relevel(dds$condition, ref = "control_NR")
#dds$group <- relevel(dds$group, ref = "control_NR")
dds$condition <- droplevels(dds$condition)
dds$group <- droplevels(dds$group)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- DESeq(dds, test="LRT", full=~Sex+group, reduced = ~ Sex)
resultsNames(dds)

control <- "coconut_NR"
treat <- "coconut_early"

# 1. Create a contrast vector called contrast_NR.
contrast <- c("group", treat, control)

# 2. Use contrast vector in the results() to extract a results table and store that to a variable called res_table_NR.
res_table <- results(dds, contrast=contrast, alpha = 0.05)

# 3. Shrink the LFC estimates using lfcShrink() and assign it back to res_table_NR.

res_table <- lfcShrink(dds, coef=paste0("group_",treat,"_vs_",control), type="apeglm")


#d <- plotCounts(dds, gene="ENSMUSG00000095687", intgroup="group", returnData=TRUE)

#ggplot(d, aes(x=group, y=count)) +  geom_point(position=position_jitter(w=0.1,h=0)) + geom_boxplot() + scale_y_log10(breaks=c(25,100,400))

vsd <- vst(dds, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Sex)
rld <- rlog(dds, blind=FALSE)
assay(rld) <- limma::removeBatchEffect(assay(rld), rld$Sex)
head(assay(vsd), 3)

plotPCA(vsd, intgroup=c("original_Diet", "duration_of_reversion_group"))
pcaData <- plotPCA(vsd, intgroup=c("group", "Sex","original_Diet","duration_of_reversion_group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=original_Diet, shape=duration_of_reversion_group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


mdata$original_Diet <- factor(mdata$original_Diet)
mdata$reversal_status <- factor(mdata$reversal_status)
mdata$duration_of_reversion_group <- factor(mdata$duration_of_reversion_group)
mdata$duration_of_reversion_months <- factor(mdata$duration_of_reversion_months)

mdata <- mdata[match(colnames(x), mdata$unique_sample_ids),]
table(rownames(mdata)==colnames(x))

#model.matrix(~current_Diet + current_Diet:duration_of_reversion_group+current_Diet:reversal_status, mdata)
#formula(~Sex + current_Diet:duration_of_reversion_group, mdata)

dds <- DESeqDataSetFromMatrix(countData = x, colData = mdata, design = formula(~original_Diet+duration_of_reversion_group))

featureData <- data.frame(gene=rownames(x))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
						
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- factor(dds$original_Diet, levels = levels(dds$original_Diet))
dds$condition <- relevel(dds$condition, ref = "control")
dds$original_Diet <- relevel(dds$original_Diet, ref = "control")
dds$condition <- droplevels(dds$condition)
dds$original_Diet <- droplevels(dds$original_Diet)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- DESeq(dds)
resultsNames(dds)

# 3. Test for the effect of reversal on the diets
design = ~ duration_of_reversion_group + original_Diet + duration_of_reversion_group:original_Diet

design(dds) <- model.matrix(~ original_Diet+ duration_of_reversion_group + Sex + duration_of_reversion_group:original_Diet,mdata)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- DESeq(dds)
resultsNames(dds)

ganot <- xanot[,c(1,2)]
temp <- as.data.frame(results(dds, name="takedown_condition_coconut_rev_vs_coconut"))
temp$Geneid <- rownames(temp)

temp_anot <- merge(temp,ganot,by="Geneid")
temp_anot <- temp_anot[order(temp_anot$padj),]

temp_anot_sig <- subset(temp_anot,temp_anot$padj<0.05 & abs(temp_anot$log2FoldChange) >1)

rownames(temp_anot_sig) <- temp_anot_sig$gene_name

row_ha = rowAnnotation(genes = anno_text(temp_anot_sig$gene_name, which = "row", gp = gpar(fontsize=10, fontface="bold")) )

Heatmap(temp_anot_sig$log2FoldChange, show_row_names=T, show_column_names=F, right_annotation = row_ha)






temp1 <- as.data.frame(results(dds, contrast = c("takedown_condition","control","coconut_rev")))
temp1$Geneid <- rownames(temp1)
temp1_anot <- merge(temp1,ganot,by="Geneid")
temp1_anot <- temp1_anot[order(temp1_anot$padj),]
res1 <- subset(temp1_anot,temp1_anot$padj<0.05 & abs(temp1_anot$log2FoldChange) >1)
rownames(res1) <- res1$gene_name
row_ha = rowAnnotation(genes = anno_text(res1$gene_name, which = "row", gp = gpar(fontsize=6, fontface="bold")) )
x11();Heatmap(res1$log2FoldChange, show_row_names=T, show_column_names=F, right_annotation = row_ha)

temp1 <- as.data.frame(results(dds, contrast = c("takedown_condition","control","coconut")))
temp1$Geneid <- rownames(temp1)
temp1_anot <- merge(temp1,ganot,by="Geneid")
temp1_anot <- temp1_anot[order(temp1_anot$padj),]
res2 <- subset(temp1_anot,temp1_anot$padj<0.05 & abs(temp1_anot$log2FoldChange) >1)
rownames(res2) <- res2$gene_name
row_ha = rowAnnotation(genes = anno_text(res2$gene_name, which = "row", gp = gpar(fontsize=6, fontface="bold")) )
x11();Heatmap(res2$log2FoldChange, show_row_names=T, show_column_names=F, right_annotation = row_ha)

