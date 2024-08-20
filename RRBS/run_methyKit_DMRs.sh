#!/bin/bash

export LC_ALL=C
module load EBModules
module load R/4.1.0-foss-2021a

echo "##### Started `date`";

Rscript /grid/beyaz/home/subhash/projects/Diet_cohorts/RRBS/scripts/methyKit_DMRs_new_final_metadata.R ;

echo "##### Completed `date`";

#qsub -N methyKit_DMRs_correctedDesignMatrix -pe threads 16 -o methyKit_DMRs_correctedDesignMatrix.out -cwd run_methyKit_DMRs.sh

