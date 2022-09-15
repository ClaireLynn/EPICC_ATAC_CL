#!/bin/bash
#
#BSUB -J "deseq[1-24]"
#BSUB -n 2
#BSUB -P DMPAVSAAO
#BSUB -W 6:00
#BSUB -o output.%J.%I
#BSUB -e errors.%J.%I
#BSUB -u clynn10@icr.ac.uk
#BSUB -R "span[hosts=1]"

# -- commands you want to execute -- #

patient=`head -n ${LSB_JOBINDEX} patients.txt | tail -1`

Rscript ATAC_regional_get_normfacs_CL3.R \
counts_recurrent/${patient}_recurrent_peaks_nf_counts.rds \
~purity+sample_type \
~purity+sample_type \
../peakcoverage/peak_coverage_data_nucleosome_free_filtered.rds \
../csaw_analysis/nf/${patient}.data.100.nf.rds \
../csaw_analysis/nf/binned.rds \
normfacs

