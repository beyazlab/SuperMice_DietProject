#!/bin/bash


export LC_ALL=C
module load EBModules

module load Python/3.7.4-GCCcore-8.3.0
module load SAMtools/1.14-GCC-10.3.0
module load BWA/0.7.17-GCC-10.2.0

genome_path="/grid/beyaz/home/subhash/dbs/mm39/indexes/bwa_mem2/";
index_name="mm39";


echo "### Started indexing with BWA mem2 `date`";
#bwameth.py index-mem2 $genome_path/${index_name}.fa ;
/grid/beyaz/home/subhash/bin/bwa-mem2/bwa-mem2 index $genome_path/${index_name}.fa ;
echo "### Completed indexing with BWA mem2 `date`";


wait

#qsub -N mm39_bwameme2_indexing -pe threads 8 -o mm39_bwameme2_indexing.out -l m_mem_free=10G -l h_rss=30G -cwd run_bwa_meth_indexing.sh 
