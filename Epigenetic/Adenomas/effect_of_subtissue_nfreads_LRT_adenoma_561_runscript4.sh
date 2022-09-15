#!/bin/bash
#
#BSUB -J "deseq"
#BSUB -n 1
#BSUB -P DMPAVSAAO
#BSUB -W 6:00
#BSUB -o output.%J.%I
#BSUB -e errors.%J.%I
#BSUB -u clynn10@icr.ac.uk
#BSUB -R "span[hosts=1]"

# -- commands you want to execute -- #


Rscript deseq_LRT_adenoma_561.R \
counts_recurrent/C561_recurrent_peaks_nf_counts.rds \
~purity+sample_type+subtissue \
~purity+sample_type \
effect_of_subtissue_nfreads_LRT \
nf
