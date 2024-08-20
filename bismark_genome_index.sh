#!/bin/bash


export LC_ALL=C
module load EBModules
module load SAMtools/1.14-GCC-10.3.0

#inpath="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/fastq";

#/grid/beyaz/home/subhash/bin/Bismark-0.23.1
#/grid/beyaz/home/subhash/bin/hisat2-2.2.1/


echo "## Bowtie2 indexing and bisulfite conversion started `date`";

/grid/beyaz/home/subhash/bin/Bismark-0.23.1/bismark_genome_preparation --parallel 8 --bowtie2 --path_to_aligner /grid/beyaz/home/subhash/bin/bowtie2-2.4.5-linux-x86_64/ /grid/beyaz/home/subhash/dbs/mm39/indexes/bowtie2_index/ ;

echo "## Bowtie2 indexing and bisulfite conversion completed `date`";


#qsub -N bismark_index -pe threads 8 -o bismark_index.out -l m_mem_free=5G -l h_rss=10G -cwd bismark_genome_index.sh

