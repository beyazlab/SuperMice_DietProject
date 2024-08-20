#export LC_ALL=C
#module load EBModules
#module load R/4.1.0-foss-2021a
#R


library(methylKit)
library(stringr)
library(vegan)
library(scales)
library(genomation)
library(corrplot)
library(BRGenomics)

inpath <- "/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/methylkit_out/"
meta_path <- "/grid/beyaz/home/subhash/projects/Diet_cohorts/metadata/"
design_matrix <- read.table(paste0(meta_path,"comparison_design_matrix_RRBS_RNAseq.txt"), header=T, sep="\t")
mdata <- read.table(paste0(meta_path,"supermice_dietcohort_complete_metadata_RRBS_RNAseq.txt"), header=T, sep="\t")
mdata$group <- paste0(mdata$original_Diet,"_",mdata$duration_of_reversion_group)


j=0
for(oname in design_matrix$output_name){
	#oname <- "coconut_early_vs_coconut_NR"
	j=j+1;
	
	inpath <- "/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/methylkit_out/"
	outpath <- paste0(inpath,"/analysis/")

	control <- subset(design_matrix$control_group,design_matrix$output_name==oname)
	subject <- subset(design_matrix$treatment_group,design_matrix$output_name==oname)
	
	cat(paste0("### Analyzing ",oname,": comparison ",j,"/",dim(design_matrix)[1] ," :",Sys.time(),"\n"))
	
	mdatax <- subset(mdata, mdata$group %in% c(control,subject))
	rownames(mdatax) <- mdatax$unique_sample_ids


	cat(paste0("### Calling methylation from diet ",oname,": ",Sys.time(),"\n"))

	file_list <- as.list(c(paste0(inpath,mdatax$unique_sample_ids,"_CpG.txt")))
	
	# Generate simple/minimal design matrix (corrected 08/20/2024)
	tmtNum <- length(subset(mdatax$group,mdatax$group==subject))
	ctlNum <- length(subset(mdatax$group,mdatax$group==control))
	design <- c(subset(mdatax[,"group"],mdatax$group==subject),subset(mdatax[,"group"],mdatax$group==control))#mdatax$group
	sid <- c(subset(mdatax[,"unique_sample_ids"],mdatax$group==subject),subset(mdatax[,"unique_sample_ids"],mdatax$group==control))#mdatax$group
	design_bin <- c(rep(1,tmtNum),rep(0,ctlNum))
	msex <- c(subset(mdatax[,"Sex"],mdatax$group==subject),subset(mdatax[,"Sex"],mdatax$group==control))#mdatax$group
	mdatay <- data.frame(sample.ids=sid,treatment=design_bin,groups=design,Sex=msex)
	
	# read the files to a methylRawList object: myobj
	myobj=methRead(file_list,
	           sample.id=as.list(c(mdatay$sample.ids)),
	           assembly="mm39",
	           treatment=c(mdatay$treatment),
	           context="CpG",
	           mincov = 10
	           )
	
	cat(paste0("### Performing differential methylation analysis for ",oname,": comparison ",j,"/",dim(design_matrix)[1] ," :",Sys.time(),"\n"))
	
	nsample <- length(myobj)
	myobjx <- filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
	newObj = normalizeCoverage(myobjx,method="median")
	
	pdf(paste0(outpath,oname,"_methylation_basicStats.pdf"), width=13, height=19)
	
	for(i in 1:nsample){
	par(mfrow=c(3,2))
	## Original calls
	print( getMethylationStats(myobj[[i]],plot=TRUE,both.strands=FALSE) )
	print( getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE) )
	## Filtered calls
	print( getMethylationStats(myobjx[[i]],plot=TRUE,both.strands=FALSE) )
	print( getCoverageStats(myobjx[[i]],plot=TRUE,both.strands=FALSE) )
	## Filtered & normalized calls
	print( getMethylationStats(newObj[[i]],plot=TRUE,both.strands=FALSE) )
	print( getCoverageStats(newObj[[i]],plot=TRUE,both.strands=FALSE) )
	}
	
	dev.off()

	methf <- unite(newObj,min.per.group=1L)
	
	pdf(paste0(outpath,oname,"_sample_correlations.pdf"), width=13, height=19)
	percmeth <- percMethylation(methf)
	percmeth[is.na(percmeth)] <- 0
	mycor <- cor(percmeth)
	colnames(mycor) <- gsub("intestine_","",colnames(mycor))
	rownames(mycor) <- gsub("intestine_","",rownames(mycor))
	corrplot(mycor,
	         title = "\n\n\n",
	         is.corr = FALSE,
	         order = "hclust",
	         col.lim = c(min(mycor),max(mycor)), 
	         col = colorRampPalette(c("white", "deepskyblue", "blue4"))(100),
	         addCoef.col = TRUE,
	         method = "shade",
	         tl.col = "black",
	         diag = TRUE)
	dev.off()
	
	pdf(paste0(outpath,oname,"_sample_clusters.pdf"), height=5)
	clusterSamples(methf, dist="correlation", method="ward.D2", plot=TRUE)
	dev.off()
	
	
	#summary(pcaData)
	
	pdf(paste0(outpath,oname,"_sample_PCAs.pdf"))
	pcaData <- PCASamples(methf, obj.return = TRUE)
	#https://github.com/fish546-2018/yaamini-virginica/blob/master/analyses/2018-10-25-MethylKit/2018-10-25-MethylKit.Rmd
	print( fig.allDataPCA <- ordiplot(pcaData, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n") )
	print( points(fig.allDataPCA, "sites", bg = alpha(c("#FA5858","#2E2E2E"),0.5), col="black", cex = 1, pch=21) )

	#Add multiple white boxes on top of the default black box to manually change the color
	print( box(col = "white") )
	print( box(col = "white") )

	print( axis(side =  1, labels = TRUE, col = "grey80", cex.axis = 1) )#Add x-axis
	print( mtext(side = 1, text = "PC 1", line = 3, cex = 1) )#Add x-axis label
	print( axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1) )#Add y-axis
	print( mtext(side = 2, text = "PC 2", line = 3, cex = 1) )#Add y-axis label
	print( text(fig.allDataPCA$sites, labels=c(gsub("intestine_","",rownames(fig.allDataPCA$sites))), cex = 0.3, offset = 0.5, pos = 2) )
	print( legend("top", 
			horiz=TRUE, inset = c(- 0.4, 0),
	       legend = unique(design), 
		   col="black",
	       fill =  alpha(c("#FA5858","#2E2E2E"),0.5), 
	       cex = 1.7, bty = "n") )#Add a legend 08with information about ambient and elevated samples
	PCASamples(methf)
	par(mfrow=c(2,2))
	print( PCASamples(methf, screeplot=TRUE) )

	dev.off()
	
	
	
	covariates=data.frame(Sex=c(mdatay$Sex))
	covariates$Sex=gsub("M","1",covariates$Sex)
	covariates$Sex=gsub("F","0",covariates$Sex)
	covariates$Sex=factor(covariates$Sex, levels = levels(factor(covariates$Sex)))
	
	myDiff=calculateDiffMeth(methf, covariates=covariates, overdispersion="MN",test="Chisq",mc.cores=16)
	
	# get hyper methylated bases
	myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
	#
	# get hypo methylated bases
	myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
	#
	#
	# get all differentially methylated bases
	myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
	
	DMRs <- as.data.frame(getData(myDiff25p[order(-abs(myDiff25p$meth.diff))]))
	hypo.DMRs <- as.data.frame(getData(myDiff25p.hypo[order(-abs(myDiff25p.hypo$meth.diff))]))
	hyper.DMRs <- as.data.frame(getData(myDiff25p.hyper[order(-abs(myDiff25p.hyper$meth.diff))]))
	#mdatay <- data.frame(sample.ids=methf@sample.ids,treatment=methf@treatment,groups=design)
	
	write.table(DMRs,paste0(outpath,oname,"_DMRs.all.txt"), quote=F, row.names=F, sep="\t")
	write.table(hyper.DMRs,paste0(outpath,oname,"_DMRs.hyper.txt"), quote=F, row.names=F, sep="\t")
	write.table(hypo.DMRs,paste0(outpath,oname,"_DMRs.hypo.txt"), quote=F, row.names=F, sep="\t")
	write.table(mdatay,paste0(outpath,oname,"_design_matrix.txt"), quote=F, row.names=F, sep="\t")
	
	gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.mm39.bed.txt", package = "methylKit"))
	

	# annotate differentially methylated CpGs with 
	# promoter/exon/intron using annotation data
	#
	# read the shores and flanking regions and name the flanks as shores 
	# and CpG islands as CpGi
	cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.mm39.bed.txt", package = "methylKit"), feature.flank.name=c("CpGi","shores"))
	#
	# convert methylDiff object to GRanges and annotate
	diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"), cpg.obj$CpGi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")
	diffCpGann.hypo=annotateWithFeatureFlank(as(myDiff25p.hypo,"GRanges"), cpg.obj$CpGi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")
	diffCpGann.hyper=annotateWithFeatureFlank(as(myDiff25p.hyper,"GRanges"), cpg.obj$CpGi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")
	
	#as.data.frame(subsetByOverlaps(as(myDiff25p,"GRanges"), c(cpg.obj$CpGi,cpg.obj$shores)))
	
	known_cpgi_DMRs=as.data.frame(subsetByOverlaps(as(myDiff25p,"GRanges"), c(cpg.obj$CpGi)))
	known_cpgi_DMRs.hypo=as.data.frame(subsetByOverlaps(as(myDiff25p.hypo,"GRanges"), c(cpg.obj$CpGi)))
	known_cpgi_DMRs.hyper=as.data.frame(subsetByOverlaps(as(myDiff25p.hyper,"GRanges"), c(cpg.obj$CpGi)))
	
	write.table(known_cpgi_DMRs,paste0(outpath,oname,"_known_cpgi_DMCs_all.txt"), quote=F, row.names=F, sep="\t")
	write.table(known_cpgi_DMRs.hyper,paste0(outpath,oname,"_known_cpgi_DMCs.hyper.txt"), quote=F, row.names=F, sep="\t")
	write.table(known_cpgi_DMRs.hypo,paste0(outpath,oname,"_known_cpgi_DMCs.hypo.txt"), quote=F, row.names=F, sep="\t")
	
	#source("/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/scripts/dataset_functions.R")
	#diffCpGannAnn=mergeGRangesData(cpg.obj,as(myDiff25p,"GRanges"), ncores=8)
	
	
	diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
	diffAnn.hypo=annotateWithGeneParts(as(myDiff25p.hypo,"GRanges"),gene.obj)
	diffAnn.hyper=annotateWithGeneParts(as(myDiff25p.hyper,"GRanges"),gene.obj)
	
	#getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

	tss_associated <- as.data.frame(getAssociationWithTSS(diffAnn))
	tss_genelist_gs <- gsub("ENS.*_","" ,unique(subset(tss_associated$feature.name,abs(tss_associated$dist.to.feature)<2000)) )
	tss_genelist_gi <- gsub("_.*","" ,unique(subset(tss_associated$feature.name,abs(tss_associated$dist.to.feature)<2000)) )
	tss_genes <- data.frame("gene_id"=tss_genelist_gi,"gene_name"=tss_genelist_gs)

	tss_associated.hypo <- as.data.frame(getAssociationWithTSS(diffAnn.hypo))
	tss_genelist_gs.hypo <- gsub("ENS.*_","" ,unique(subset(tss_associated.hypo$feature.name,abs(tss_associated.hypo$dist.to.feature)<2000)) )
	tss_genelist_gi.hypo <- gsub("_.*","" ,unique(subset(tss_associated.hypo$feature.name,abs(tss_associated.hypo$dist.to.feature)<2000)) )
	tss_genes.hypo <- data.frame("gene_id"=tss_genelist_gi.hypo,"gene_name"=tss_genelist_gs.hypo)
	
	tss_associated.hyper <- as.data.frame(getAssociationWithTSS(diffAnn.hyper))
	tss_genelist_gs.hyper <- gsub("ENS.*_","" ,unique(subset(tss_associated.hyper$feature.name,abs(tss_associated.hyper$dist.to.feature)<2000)) )
	tss_genelist_gi.hyper <- gsub("_.*","" ,unique(subset(tss_associated.hyper$feature.name,abs(tss_associated.hyper$dist.to.feature)<2000)) )
	tss_genes.hyper <- data.frame("gene_id"=tss_genelist_gi.hyper,"gene_name"=tss_genelist_gs.hyper)
	
	
	write.table(tss_genes,paste0(outpath,oname,"_TSS_methylated_genelist.all.txt"), quote=F, row.names=F, sep="\t")
	write.table(tss_genes.hyper,paste0(outpath,oname,"_TSS_methylated_genelist.hyper.txt"), quote=F, row.names=F, sep="\t")
	write.table(tss_genes.hypo,paste0(outpath,oname,"_TSS_methylated_genelist.hypo.txt"), quote=F, row.names=F, sep="\t")
	
	write.table(unique(tss_genes$gene_name),paste0(outpath,"functional_analysis/input/",oname,"_TSS_methylated_genes.all.list"), quote=F, row.names=F, sep="\n", col.names=F)
	write.table(unique(tss_genes.hyper$gene_name),paste0(outpath,"functional_analysis/input/",oname,"_TSS_methylated_genes.hyper.list"), quote=F, row.names=F, sep="\n", col.names=F)
	write.table(unique(tss_genes.hypo$gene_name),paste0(outpath,"functional_analysis/input/",oname,"_TSS_methylated_genes.hypo.list"), quote=F, row.names=F, sep="\n", col.names=F)
	
	write.table(tss_associated,paste0(outpath,oname,"_TSS_methylated_genes_distance.all.txt"), quote=F, row.names=F, sep="\t")
	write.table(tss_associated.hyper,paste0(outpath,oname,"_TSS_methylated_genes_distance.hyper.txt"), quote=F, row.names=F, sep="\t")
	write.table(tss_associated.hypo,paste0(outpath,oname,"_TSS_methylated_genes_distance.hypo.txt"), quote=F, row.names=F, sep="\t")
	

	pdf(paste0(outpath,oname,"_diffMeth_stats.pdf"), width=9, height=12)
	par(mfrow=c(2,2))
	plotTargetAnnotation(diffAnn,precedence=TRUE, main=paste0("DMRs all: ",oname))
	plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),  main=paste0("DMRs all: ",oname))
	barplot(getFeatsWithTargetsStats(diffAnn,percentage=TRUE), col="#424242")
	par(mfrow=c(2,2))
	plotTargetAnnotation(diffAnn.hyper,precedence=TRUE, main=paste0("DMRs hyper: ",oname))
	plotTargetAnnotation(diffCpGann.hyper,col=c("green","gray","white"),  main=paste0("DMRs hyper: ",oname))
	barplot(getFeatsWithTargetsStats(diffAnn.hyper,percentage=TRUE), col="#424242")
	par(mfrow=c(2,2))
	plotTargetAnnotation(diffAnn.hypo,precedence=TRUE, main=paste0("DMRs hypo: ",oname))
	plotTargetAnnotation(diffCpGann.hypo,col=c("green","gray","white"),  main=paste0("DMRs hypo: ",oname))
	barplot(getFeatsWithTargetsStats(diffAnn.hypo,percentage=TRUE), col="#424242")
	
	dev.off()
	
	saveRDS(methf, paste0(outpath,oname,"_filtered_methObj.rds"))

	cat(paste0("### Completed differential methylation analysis for ",oname,": ",Sys.time(),"\n"))

	
	
}


