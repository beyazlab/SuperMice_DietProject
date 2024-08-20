library(methylKit)
library(stringr)
library(vegan)
library(scales)
library(genomation)


diets <- c("coco","fish","keto","lard","milk","olive","palm")
ctl <- "cont"

inpath <- "/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/bismark_results"
outpath <- "/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/methylkit_out/"

for(tmt in diets){
	
	cat(paste0("### Calling methylation from diet ",tmt,": ",Sys.time(),"\n"))
	
	file_list_c <- as.list(list.files(path = inpath, recursive = TRUE, pattern = paste0("intestine_",ctl,"_.*[0-9].bam$"), full.names = TRUE))
	file_list_t <- as.list(list.files(path = inpath, recursive = TRUE, pattern = paste0("intestine_",tmt,"_.*[0-9].bam$"), full.names = TRUE))
	intestine_file_list <- c(file_list_c,file_list_t)
	sample_names <- list.files(path = inpath, recursive = TRUE, pattern = paste0("[0-9].bam$"), full.names = FALSE)
	intestine_sample_names <- str_extract(sample_names, paste(c(paste0("intestine.*",ctl,".*"),paste0("intestine.*",tmt,".*")), collapse ="|"))
	intestine_sample_names <- intestine_sample_names[!is.na(intestine_sample_names)]
	intestine_sample_names <- gsub(".bam","",intestine_sample_names)


	temp1 <- gsub("_[0-9]_[0-9][0-9]","",intestine_sample_names)
	temp2 <- gsub("_[0-9]_[0-9]","",temp1)
	temp3 <- gsub("_[0-9]_[0-9].[0-9]","",temp2)
	temp4 <- gsub(".[0-9]","",temp3)
	temp5 <- gsub("_[0-9]_[0-9][0-9]","",temp4)
	design <- gsub("intestine_","",temp5)

	design_bin <- design
	design_bin <- gsub(ctl,"0",gsub(tmt,"1", design_bin))
	design_bin <- as.numeric(design_bin)


	my.methRaw <- processBismarkAln(location=intestine_file_list, sample.id=as.list(intestine_sample_names), "mm39", save.folder = "/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/methylkit_out",
	  save.context = c("CpG"), read.context = "CpG", nolap = FALSE,
	  mincov = 10, minqual = 20, phred64 = FALSE, treatment=c(design_bin),
	  save.db = FALSE)
	
  	cat(paste0("### Methylation calling finished for diet ",tmt,": ",Sys.time(),"\n"))
  	cat(paste0("### Saving methylation calls for diet ",tmt,": ",Sys.time()))
	

	save.image(paste0(outpath,tmt,"_vs_",ctl,"_methRaw_object.Rdata"))
	
	cat(paste0("### Finished saving methylation calls for diet ",tmt,": ",Sys.time(),"\n"))
	

}

