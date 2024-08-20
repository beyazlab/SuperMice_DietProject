#!/bin/bash

#### Submit multipe jobs
#find . -name '*_trim.sh' | xargs -I% qsub -pe threads 8 -l m_mem_free=5G -l h_rss=10G -cwd %
#find . -name '*_align.sh' | xargs -I% qsub -pe threads 8 -l m_mem_free=5G -l h_rss=10G -cwd %
#find . -name '*_bismark.sh' | xargs -I% qsub -pe threads 16 -l m_mem_free=5G -l h_rss=10G -cwd %
#find . -name 'intestine_cont_*_bismark.sh' | xargs -I% qsub -pe threads 16 -cwd %

task=$1;
threads=$2;
SCRIPT_PATH="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/scripts";

index_name="mm39";
genome_path="/grid/beyaz/home/subhash/dbs/mm39/indexes/bwa_mem2";

SAMPLES=(`cat $SCRIPT_PATH/sample_metadata.txt | grep -v "^#"| cut -f1,2,3,4 | awk '!x[$0]++' | sed 's/\t/;/g'`); 


for i in ${SAMPLES[@]}

do

	IFS=';' read -ra arr <<< "$i"; 
	path=${arr[0]}
	dir=${arr[1]};
	sname=${arr[2]};
	oname=${arr[3]};

	if [ "$task" == "trim" ]; then
		inpath=$path;
		outpath="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/fastq";
		qcpath="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/QC/clean/";
		echo "#!/bin/bash" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> bscripts/${arr[3]}_${task}.sh
		echo "export LC_ALL=C" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "module load EBModules" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "alias fastqc='/grid/beyaz/home/subhash/bin/FastQC/fastqc';" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "echo \"### Started trimming $sname with Trimgalore \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		echo "" >> bscripts/${arr[3]}_${task}.sh
		echo "/grid/beyaz/home/subhash/bin/TrimGalore-0.6.6/trim_galore --output_dir ${outpath} --paired --rrbs --cores $threads --path_to_cutadapt /grid/beyaz/home/subhash/.local/bin/cutadapt --basename ${oname} --fastqc_args \"--threads $threads --outdir ${qcpath}\" $inpath/$dir/${sname}_R1_001.fastq.gz $inpath/$dir/${sname}_R2_001.fastq.gz ;" >> bscripts/${arr[3]}_${task}.sh 
		echo "echo \"### Completed trimming $sname with Trimgalore \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		echo "wait" >> bscripts/${arr[3]}_${task}.sh
		echo "" >> bscripts/${arr[3]}_${task}.sh
	fi
	if [ "$task" == "align" ]; then
		inpath="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/fastq";
		outpath="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/aligned";
		echo "#!/bin/bash" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> bscripts/${arr[3]}_${task}.sh
		echo "export LC_ALL=C" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "module load EBModules" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "module load Python/3.7.4-GCCcore-8.3.0" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "module load SAMtools/1.14-GCC-10.3.0" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "echo \"### Started aligning $sname with BWA mem2 \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		echo "" >> bscripts/${arr[3]}_${task}.sh
		echo "/grid/beyaz/home/subhash/bin/bwa-mem2/bwa-mem2 mem -L 25 -pCM -t $threads $genome_path/${index_name}.fa '<python /grid/beyaz/home/subhash/.local/bin/bwameth.py c2t $inpath/${oname}_val_1.fq.gz $inpath/${oname}_val_2.fq.gz' | samtools view -bS - | samtools sort --write-index -@ $threads - -o $outpath/${oname}.bam ;" >> bscripts/${arr[3]}_${task}.sh
		echo "echo \"## Indexing BAM \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		echo "samtools index $outpath/${oname}.bam $outpath/${oname}.bam.bai ;" >> bscripts/${arr[3]}_${task}.sh
		echo "echo \"### Completed aligning $sname with BWA mem2 \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		echo "wait" >> bscripts/${arr[3]}_${task}.sh
		echo "" >> bscripts/${arr[3]}_${task}.sh
	fi
	if [ "$task" == "bismark" ]; then
		inpath="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/fastq";
		outpath="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/bismark_results";
		echo "#!/bin/bash" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> bscripts/${arr[3]}_${task}.sh
		echo "export LC_ALL=C" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "module load EBModules" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		#echo "module load Python/3.7.4-GCCcore-8.3.0" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "module load SAMtools/1.14-GCC-10.3.0" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "" >> $SCRIPT_PATH/bscripts/${arr[3]}_${task}.sh
		echo "echo \"### Started processing of $sname with Bismark Bowtie2 \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		echo "" >> bscripts/${arr[3]}_${task}.sh
		echo "/grid/beyaz/home/subhash/bin/Bismark-0.23.1/bismark --nucleotide_coverage --output_dir ${outpath} --rg_id ${oname} --rg_sample ${oname} --parallel $threads --path_to_bowtie2 /grid/beyaz/home/subhash/bin/bowtie2-2.4.5-linux-x86_64/ --samtools_path /grid/it/data/elzar/easybuild/software/SAMtools/1.14-GCC-10.3.0/bin/ --genome /grid/beyaz/home/subhash/dbs/mm39/indexes/bowtie2_index/ -1 $inpath/${oname}_val_1.fq.gz -2 $inpath/${oname}_val_2.fq.gz ;" >> bscripts/${arr[3]}_${task}.sh
		echo "echo \"## Sorting & Indexing BAM \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		#echo "samtools view -bS $outpath/${oname}_val_1_bismark_bt2_pe.bam | samtools sort -n --write-index -@ $threads - -o $outpath/${oname}_rSorted.bam ;" >> bscripts/${arr[3]}_${task}.sh
		#echo "samtools index $outpath/${oname}_rSorted.bam $outpath/${oname}_rSorted.bam.bai ;" >> bscripts/${arr[3]}_${task}.sh

		echo "samtools sort -@ $threads $outpath/${oname}_val_1_bismark_bt2_pe.bam -o $outpath/${oname}.bam ;" >> bscripts/${arr[3]}_${task}.sh
		echo "samtools index $outpath/${oname}.bam $outpath/${oname}.bam.bai ;" >> bscripts/${arr[3]}_${task}.sh
		echo "rm $outpath/${oname}_val_1_bismark_bt2_pe.bam;" >> bscripts/${arr[3]}_${task}.sh
		#echo "echo \"## Methylation extraction started \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		#echo "/grid/beyaz/home/subhash/bin/Bismark-0.23.1/bismark_methylation_extractor --ucsc --cytosine_report --bedGraph --parallel $threads --samtools_path /grid/it/data/elzar/easybuild/software/SAMtools/1.14-GCC-10.3.0/bin/ --output ${outpath} --paired-end --genome_folder /grid/beyaz/home/subhash/dbs/mm39/indexes/bowtie2_index/ $outpath/${oname}_rSorted.bam ;" >> bscripts/${arr[3]}_${task}.sh
		#echo "echo \"## Generating reports on methylation extraction \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		#echo "/grid/beyaz/home/subhash/bin/Bismark-0.23.1/bismark2report --output ${oname}.html --dir $outpath --alignment_report $outpath/${oname}_val_1_bismark_bt2_PE_report.txt --splitting_report $outpath/${oname}_splitting_report.txt --mbias_report $outpath/${oname}.M-bias.txt --nucleotide_report $outpath/${oname}_val_1_bismark_bt2_pe.nucleotide_stats.txt;" >> bscripts/${arr[3]}_${task}.sh
		#echo "echo \"## Generating summary on methylation extraction \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		#echo "mv $outpath/${oname}_val_1_bismark_bt2_PE_report.txt $outpath/${oname}_PE_report.txt;" >> bscripts/${arr[3]}_${task}.sh
		#echo "/grid/beyaz/home/subhash/bin/Bismark-0.23.1/bismark2summary --basename $outpath/${oname}_bismark2summary $outpath/${oname}_val_1_bismark_bt2_pe.bam ;" >> bscripts/${arr[3]}_${task}.sh
		echo "echo \"### Completed processing of $sname with Bismark Bowtie2  \`date\`\";" >> bscripts/${arr[3]}_${task}.sh
		echo "wait" >> bscripts/${arr[3]}_${task}.sh
		echo "" >> bscripts/${arr[3]}_${task}.sh
	fi
	
	#qsub -pe threads 16 -l m_mem_free=5G -l h_rss=10G -cwd intestine_coco_1_1_bismark.sh
done





