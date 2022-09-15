#!/bin/bash
#
#BSUB -J "counts[24]"
#BSUB -n 2
#BSUB -P DMPAVSAAO
#BSUB -W 4:00
#BSUB -o output.%J.%I
#BSUB -e errors.%J.%I
#BSUB -u clynn10@icr.ac.uk
#BSUB -R "span[hosts=1]"

# -- commands you want to execute -- #

patient=`head -n ${LSB_JOBINDEX} patients.txt | tail -1`

Rscript ATAC_regional_get_recurrent_peaks_readcounts_CL2.R -p $patient -r recurrent_peaks.rds 
