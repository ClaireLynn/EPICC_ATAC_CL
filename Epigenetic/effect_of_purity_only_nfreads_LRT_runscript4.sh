#!/bin/bash
#
#BSUB -J "deseq[1-24]"
#BSUB -n 1
#BSUB -P DMPAVSAAO
#BSUB -W 6:00
#BSUB -o output.%J.%I
#BSUB -e errors.%J.%I
#BSUB -u clynn10@icr.ac.uk
#BSUB -R "span[hosts=1]"

# -- commands you want to execute -- #

patient=`head -n ${LSB_JOBINDEX} patients.txt | tail -1`

Rscript ATAC_regional_recurrent_peaks_deseq_LRT_CL4.R \
counts_recurrent/${patient}_recurrent_peaks_nf_counts.rds \
~sample_type+purity \
~sample_type \
effect_of_purity_only_nfreads_LRT \
nf
