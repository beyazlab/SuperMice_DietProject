#!/bin/bash


export LC_ALL=C
module load EBModules
module load R/4.1.0-foss-2021a

#inpath="/grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/fastq";

#/grid/beyaz/home/subhash/bin/Bismark-0.23.1
#/grid/beyaz/home/subhash/bin/hisat2-2.2.1/


echo "#### Generating methyKit objects started `date`";

Rscript /grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/scripts/make_methylKit_objects.R ;

echo "#### Generating methyKit objects completed `date`";


#qsub -N prepare_methylkit -pe threads 4 -o prepare_methylkit.out -cwd make_methylKit_objects.sh

