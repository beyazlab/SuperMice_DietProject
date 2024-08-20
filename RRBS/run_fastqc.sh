#!/bin/bash
export LC_ALL=C
#module load EBModules-LegacyBNB
#module load FastQC/0.11.8-Java-1.8

echo "#### Started `date`";

/grid/beyaz/home/subhash/bin/FastQC/fastqc --outdir /grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/QC/ --threads 16 /mnt/grid/beyaz/hpc/home/data/meydan/Diet_RRBS/*/*.fastq.gz ;
echo "## Batch 2 `date`";
/grid/beyaz/home/subhash/bin/FastQC/fastqc --outdir /grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/QC/ --threads 6 /mnt/grid/beyaz/hpc/home/data/meydan/Diet_RRBS_2/*/*.fastq.gz ;

echo "#### Completed `date`";
wait

#qsub -N fastqc_check -pe threads 16 -o fastqc_check.out -l m_mem_free=5G -l h_rss=10G -cwd run_fastqc.sh