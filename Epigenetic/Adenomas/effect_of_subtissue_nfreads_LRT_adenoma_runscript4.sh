#!/bin/bash
#
#BSUB -J "deseq[1-6]"
#BSUB -n 2
#BSUB -P DMPAVSAAO
#BSUB -W 6:00
#BSUB -o output.%J.%I
#BSUB -e errors.%J.%I
#BSUB -u clynn10@icr.ac.uk
#BSUB -R "span[hosts=1]"

# -- commands you want to execute -- #

patient=`head -n ${LSB_JOBINDEX} adenoma_patients.txt | tail -1`

Rscript deseq_LRT_adenoma.R \
counts_recurrent/${patient}_recurrent_peaks_nf_counts.rds \
~purity+sample_type+subtissue \
~purity+sample_type \
effect_of_subtissue_nfreads_LRT \
nf
