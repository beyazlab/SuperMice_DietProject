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
mdata$group <- paste0(mdata$original_Diet,"_",mdata$duration_of_reversion_group)
design_matrix <- read.table(paste0(meta_path,"comparison_design_matrix_RRBS_RNAseq.txt"), header=T, sep="\t")

#control <- "control_NR"
#subject <- "coconut_early"
#oname <- paste0(subject,"_vs_",control)
lfc <- 1
fdr <- 0.05
i=0

for(oname in design_matrix$output_name){
	i=i+1;
x <- read.table(paste0(inpath,"merged_gene_counts_matchingRRBS_clean.txt"), header=T, sep="\t", row.names="Geneid")
control <- subset(design_matrix$control_group,design_matrix$output_name==oname)
subject <- subset(design_matrix$treatment_group,design_matrix$output_name==oname)

cat(paste0("### Analyzing ",control,"_vs_",subject,": comparison ",i,"/",dim(design_matrix)[1] ," :",Sys.time(),"\n"))

mdatax <- subset(mdata, mdata$group %in% c(control,subject))
rownames(mdatax) <- mdatax$unique_sample_ids
x <- x[,mdatax$unique_sample_ids]

#### Common test

dds <- DESeqDataSetFromMatrix(countData = x, colData = mdatax, design =~Sex+group)
featureData <- data.frame(gene=rownames(x))
mcols(dds) <- DataFrame(mcols(dds), featureData)
#mcols(dds)
#keep <- rowSums(counts(dds)) >= 10
keep <- rowSums(counts(dds) >= 10) >= dim(mdatax)[1]/2
dds <- dds[keep,]
dds$condition <- factor(dds$group, levels = levels(dds$group))
dds$condition <- relevel(dds$condition, ref = control)
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- DESeq(dds, test="LRT", full=~Sex+group, reduced = ~ Sex)
resultsNames(dds)
#norm_counts <- counts(dds, normalized=TRUE)
norm_counts <- limma::removeBatchEffect(log2(counts(dds, normalized=TRUE)+1), dds$Sex)

norm_counts_anotd <- merge(xanot[,c("Geneid","gene_name")],norm_counts,by.x="Geneid",by.y="row.names")
write.table(norm_counts_anotd, paste0(outpath,oname,"_normCounts.txt"), sep="\t", quote=F, row.names=F)
rownames(norm_counts_anotd) <- norm_counts_anotd$gene_id

vsd <- vst(dds, blind=FALSE)

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Sex)

pdf(paste0(plot_path,oname,"_PCA_plots.pdf"))
p <- pca(assay(vsd), metadata = mdatax, removeVar = 0.1)
print( screeplot(p, axisLabSize = 18, titleLabSize = 22) )
print( biplot(p, showLoadings = TRUE, lab = NULL) )
#pairsplot(p)
horn <- parallelPCA(assay(vsd))
#horn$n
elbow <- findElbowPoint(p$variance)
#elbow
print( screeplot(p,components = getComponents(p, 1:20),vline = c(horn$n, elbow)) + geom_label(aes(x = horn$n + 1, y = 50,label = 'Horn\'s', vjust = -1, size = 8)) + geom_label(aes(x = elbow + 1, y = 50,label = 'Elbow method', vjust = -1, size = 8)) )
#which(cumsum(p$variance) > 80)[1]
print( biplot(p, lab = paste0(p$metadata$group), colby = 'original_Diet', shape='Sex', hline = 0, vline = 0, legendPosition = 'right') )
print( biplot(p, lab = paste0(p$metadata$group), colby = 'cage', shape='Sex', hline = 0, vline = 0, legendPosition = 'right') )

print( biplot(p, lab = paste0(p$metadata$group), colby = 'duration_of_reversion_group', hline = 0, vline = 0, legendPosition = 'right') )
print( biplot(p, lab = paste0(p$metadata$group), colby = 'Sex', hline = 0, vline = 0, legendPosition = 'right') )

dev.off()


pdf(paste0(plot_path,oname,"_Overall_transformation_plots.pdf"));
ntd <- normTransform(dds)

print( meanSdPlot(assay(ntd)) ) ## Non-transformed data
print( meanSdPlot(assay(vsd)) ) ## Transformed data
dev.off()

res <- as.data.frame(results(dds, contrast=c("group", subject, control), independentFiltering=TRUE, alpha=fdr, pAdjustMethod="BH", parallel=F))
res <- merge(res, norm_counts_anotd, by.x="row.names",by.y="Geneid")
colnames(res)[1] <- "gene_id"
res_sig <- subset(res, abs(res$log2FoldChange)>lfc & res$padj<fdr)
res_sig <- res_sig[order(res_sig$padj),]
up <- dim(subset(res_sig,res_sig$log2FoldChange>lfc ))[1]
down <- dim(subset(res_sig,res_sig$log2FoldChange<lfc ))[1]
DEup <- subset(res_sig,res_sig$log2FoldChange>lfc )
DEdown <- subset(res_sig,res_sig$log2FoldChange<lfc )
write.table(DEup$gene_name,paste0(frich_path,"input/",oname,"_DE_up.list"), col.names=F, row.names=F, sep="\n", quote=F)
write.table(DEdown$gene_name,paste0(frich_path,"input/",oname,"_DE_down.list"), col.names=F, row.names=F, sep="\n", quote=F)
write.table(c(DEup$gene_name,DEdown$gene_name),paste0(frich_path,"input/",oname,"_DE_all.list"), col.names=F, row.names=F, sep="\n", quote=F)
write.table(res_sig, paste0(outpath,oname,"_DE_sig_genes.txt"), quote=F, sep="\t", row.names=F)
write.table(res, paste0(outpath,oname,"_DE_allStats.txt"), quote=F, sep="\t", row.names=F)


if(dim(res_sig)[1]>1){

pdf(paste0(plot_path,oname,"_DE_genes_heatmap.pdf"))
temp <- res_sig
snames <- mdatax$unique_sample_ids
expMatx <- temp[,c("gene_name",snames)] #,"baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","Geneid","gene_name"
temp$GeneSymbol <- make.names(temp$gene_name, unique = TRUE)
rownames(expMatx) <- temp$GeneSymbol
expMatx <- expMatx[,2:ncol(expMatx)]
mat <- t(apply(expMatx, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))
mat1 <- mat[,grepl("WT", colnames(mat))]
mat2 <- mat[,!(grepl("WT", colnames(mat)))]
hclustfunc <- function(x, method = "complete", dmeth = "euclidean") {    
    hclust(dist(x, method = dmeth), method = method)
}
col_fun = colorRamp2(c(min(mat), 0, max(mat)), c("#0404B4", "white", "#DF013A"))
col_fun1 = colorRamp2(c(min(mat1), 0, max(mat1)), c("#0404B4", "white", "#DF013A"))
col_fun2 = colorRamp2(c(min(mat2), 0, max(mat2)), c("#0404B4", "white", "#DF013A"))
km = kmeans(mat, centers = 2)$cluster
row_order = hclustfunc(mat)$order

ht1 <- Heatmap(mat1, col=col_fun1, row_order = row_order, cluster_rows = TRUE, split = km, name=paste0(control), column_title=paste0(oname," DEGs"), column_names_gp = gpar(fontsize = 5)) + Heatmap(mat2, col=col_fun2, name=paste0(subject), row_names_gp = gpar(fontsize = 2), column_names_gp = gpar(fontsize = 5)) 
ht2 <- Heatmap(mat, col=col_fun, cluster_rows = TRUE, km=2, name=paste0("Expr"), column_title=paste0(oname," DEGs"), column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 2))
ht3 <- Heatmap(mat, col=col_fun, column_split=mdatax$group, cluster_rows = TRUE, km=2, name=paste0("Expr"), column_title=paste0(oname," DEGs"), column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 2))

draw(ht1)
draw(ht2)
draw(ht3)


dev.off()


pdf(paste0(plot_path,oname,"_DE_genes_volcano_all.pdf"))
temp1 <- res
#subset(NR_vs_CC_sig, NR_vs_CC_sig$gene_biotype == "protein_coding")

temp1$GeneSymbol <- make.names(temp1$gene_name, unique = TRUE)
rownames(temp1) <- temp1$GeneSymbol
de <- temp1
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < fdr ] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < fdr ] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- rownames(de[de$diffexpressed != "NO",])

p1 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
        geom_point(size = 2) + 
        theme(axis.text = element_text(size = 6, color="black"), axis.title = element_text(size = 6), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "#2E2E2E", fill=NA, size=0.5)) +
        geom_text_repel(size=2, fontface="bold", box.padding = 0.2, segment.size = 0.2, min.segment.length = 0.02, segment.color = '#585858', segment.curvature = 0, segment.ncp = 1, segment.angle = 20, arrow=arrow(angle = 20, length = unit(0.01, "inches"), ends = "last", type = "open")) +
        scale_color_manual(values=c("#5F04B4","#424242", "#B4045F")) +
        geom_vline(xintercept=c(-lfc, lfc), col="#585858", linetype="longdash", size=0.3) +
        geom_hline(yintercept=-log10(fdr), col="#585858", linetype="longdash", size=0.3) + ggtitle(paste0(oname)) + xlab(paste0("log2(",subject,"/",control,",)")) + ylab("-log10(FDR)")
print(p1)
dev.off()



## gene plots
if(dim(res_sig)[1]>5){
genes <- res_sig$gene_name[1:6]
gene_mat <- subset(norm_counts_anotd, norm_counts_anotd$gene_name %in% c(genes))
gene_mat <- gene_mat[,c("gene_name",mdatax$unique_sample_ids)]
rownames(mdatax) <- mdatax$unique_sample_ids
rownames(gene_mat) <- gene_mat$gene_name
gene_mat <- gene_mat[,2:length(colnames(gene_mat))]
gene_matx <- t(gene_mat)
gene_matx <- merge(gene_matx,mdatax,by="row.names")
gene_matx <- as.data.frame(gene_matx)
gene_matx$group <- gsub(subject,paste0("z",subject),gene_matx$group)

pdf(paste0(outpath,"plots/",oname,"_top_genes_expr.pdf"))

p1 <- ggplot(gene_matx,aes(gene_matx[,genes[1]],group, fill=group))+ylab("")+xlab(paste0(genes[1]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))
p2 <- ggplot(gene_matx,aes(gene_matx[,genes[2]],group, fill=group))+ylab("")+xlab(paste0(genes[2]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))
p3 <- ggplot(gene_matx,aes(gene_matx[,genes[3]],group, fill=group))+ylab("")+xlab(paste0(genes[3]," expression"))+geom_boxplot()+coord_flip()+theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))

p4 <- ggplot(gene_matx,aes(log2(gene_matx[,genes[1]]),group, fill=group))+ylab("")+xlab(paste0(genes[1]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))
p5 <- ggplot(gene_matx,aes(log2(gene_matx[,genes[2]]),group, fill=group))+ylab("")+xlab(paste0(genes[2]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))
p6 <- ggplot(gene_matx,aes(log2(gene_matx[,genes[3]]),group, fill=group))+ylab("")+xlab(paste0(genes[3]," expression"))+geom_boxplot()+coord_flip()+theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))

print( p1+p2+p3+p4+p5+p6+plot_layout(ncol=3,nrow=2)+ plot_annotation(title = paste0(oname),theme = theme(plot.title = element_text(hjust = 0.5))) )

p1 <- ggplot(gene_matx,aes(gene_matx[,genes[4]],group, fill=group))+ylab("")+xlab(paste0(genes[4]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))
p2 <- ggplot(gene_matx,aes(gene_matx[,genes[5]],group, fill=group))+ylab("")+xlab(paste0(genes[5]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))
p3 <- ggplot(gene_matx,aes(gene_matx[,genes[6]],group, fill=group))+ylab("")+xlab(paste0(genes[6]," expression"))+geom_boxplot()+coord_flip()+theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))

p4 <- ggplot(gene_matx,aes(log2(gene_matx[,genes[4]]),group, fill=group))+ylab("")+xlab(paste0(genes[4]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))
p5 <- ggplot(gene_matx,aes(log2(gene_matx[,genes[5]]),group, fill=group))+ylab("")+xlab(paste0(genes[5]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))
p6 <- ggplot(gene_matx,aes(log2(gene_matx[,genes[6]]),group, fill=group))+ylab("")+xlab(paste0(genes[6]," expression"))+geom_boxplot()+coord_flip()+theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))


print( p1+p2+p3+p4+p5+p6+plot_layout(ncol=3,nrow=2)+ plot_annotation(title = paste0(oname),theme = theme(plot.title = element_text(hjust = 0.5))) )

dev.off()
}
if(dim(res_sig)[1]<5){

	genes <- res_sig$gene_name[1:2]
	gene_mat <- subset(norm_counts_anotd, norm_counts_anotd$gene_name %in% c(genes))
	gene_mat <- gene_mat[,c("gene_name",mdatax$unique_sample_ids)]
	rownames(mdatax) <- mdatax$unique_sample_ids
	rownames(gene_mat) <- gene_mat$gene_name
	gene_mat <- gene_mat[,2:length(colnames(gene_mat))]
	gene_matx <- t(gene_mat)
	gene_matx <- merge(gene_matx,mdatax,by="row.names")
	gene_matx <- as.data.frame(gene_matx)
	gene_matx$group <- gsub(subject,paste0("z",subject),gene_matx$group)

	pdf(paste0(outpath,"plots/",oname,"_top_genes_expr.pdf"))

	p1 <- ggplot(gene_matx,aes(gene_matx[,genes[1]],group, fill=group))+ylab("")+xlab(paste0(genes[1]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))
	p2 <- ggplot(gene_matx,aes(gene_matx[,genes[2]],group, fill=group))+ylab("")+xlab(paste0(genes[2]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))

	p3 <- ggplot(gene_matx,aes(log2(gene_matx[,genes[1]]),group, fill=group))+ylab("")+xlab(paste0(genes[1]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))
	p4 <- ggplot(gene_matx,aes(log2(gene_matx[,genes[2]]),group, fill=group))+ylab("")+xlab(paste0(genes[2]," expression"))+geom_boxplot()+coord_flip()+theme(legend.position = "none", axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#585858","#04B45F"))

	print( p1+p2+p3+p4+plot_layout(ncol=2,nrow=2)+ plot_annotation(title = paste0(oname),theme = theme(plot.title = element_text(hjust = 0.5))) )

	dev.off()
	
}

ghighlight <- res_sig$gene_name[1:round(length(res_sig$gene_name)*0.40)]
DEg<- res_sig
rownames(DEg) <- DEg$gene_name

matx <- DEg[,9:length(colnames(DEg))]

mat <- t(apply(matx, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))

col_fun = colorRamp2(c(min(mat), 0, max(mat)), c("#08088A", "white", "#B40431")) #c("#0B3861", "white", "#B40431")

row_ha = rowAnnotation(foo = anno_mark(at = c(which(rownames(mat) %in% ghighlight)), labels = subset(rownames(mat),rownames(mat) %in% ghighlight), labels_gp = gpar(fontsize=8,fontface = "italic")), annotation_name_gp= gpar(fontsize = 2) )

ht1<- Heatmap(mat,column_title=paste0(oname), show_column_dend = FALSE, heatmap_legend_param = list(legend_direction = "vertical"), show_row_dend = FALSE, row_names_rot = 90, column_names_rot = -90, right_annotation = row_ha, col=col_fun, row_title_gp=gpar(fontsize=1), column_title_gp=gpar(fontsize=8), column_split=mdatax$group, show_row_names=F, cluster_rows = TRUE, km=2, name=paste0("Expr"), column_names_gp = gpar(fontsize = 8), row_title=NULL)

pdf(paste0(plot_path,oname,"_DE_heatmap_top_genes_highlighted.pdf"))

gb = grid.grabExpr(draw(ht1, heatmap_legend_side = "left", annotation_legend_side = "left"))
gb
is.grob(gb)
pushViewport(viewport(width = 0.4, height = 1,angle = 90))
grid.draw(gb)

dev.off()
}

}










