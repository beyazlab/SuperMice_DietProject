#!/bin/bash


export LC_ALL=C
module load EBModules

sname=$1;
#sname="Control_1-1-F_SF_Ileum_7_7_20_S1_L002";
dir=$2 #"EC-CM-6649_1";
oname=$3 #"Control_1-1";

inpath="/mnt/grid/beyaz/hpc/home/data/meydan/Diet_RRBS";
qcpath="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/QC/clean/";
outpath="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/fastq/"

alias fastqc='/grid/beyaz/home/subhash/bin/FastQC/fastqc';

echo "#### Started `date`";

/grid/beyaz/home/subhash/bin/TrimGalore-0.6.6/trim_galore --output_dir ${outpath} --paired --rrbs --cores 8 --path_to_cutadapt /grid/beyaz/home/subhash/.local/bin/cutadapt --basename ${oname} --fastqc_args "--threads 2 --outdir ${qcpath}" $inpath/$dir/${sname}_R1_001.fastq.gz $inpath/$dir/${sname}_R2_001.fastq.gz ;

echo "#### Completed `date`";
wait

#qsub -N trim_galore_test -pe threads 8 -o trim_galore_test.out -l m_mem_free=5G -l h_rss=10G -cwd run_trim_galore.sh
