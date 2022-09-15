# This script gets the coverage of enhancers for specified samples (WT/mutant),
# performs normalisation per sample and performs a T-Test

library(rtracklayer)

# read peak counts and matching peaks for each patient
counts <- readRDS("atacseq_counts_per_peak.rds")
peak_matches_df<-readRDS("peak_matches_df.rds")
colnames(counts) <- gsub("EPICC_|_C1$","",colnames(counts))

# get eqtl results and matching WGS-ATAC samples
eqtl<- readRDS("alternative_eqtls_plus_HWrec.newclo.rds")
matches <- readRDS("dna_sample_data_with_atac.rds")
samples <-unique(unlist(strsplit(eqtl$MutAllSam,";")))
patients <- unique(gsub("_.*","",samples))

# get mutation positions
eqtl$chr <- gsub(":.*","",eqtl$ID)
eqtl$start <-  gsub("_.*","",gsub(".*:","",eqtl$ID))
eqtl$end <- eqtl$start
muts <- makeGRangesFromDataFrame(eqtl,
                                 keep.extra.columns=T,
                                 ignore.strand=FALSE,
                                 seqinfo=NULL,
                                 starts.in.df.are.0based=FALSE)


# do T-test between mutant samples and wt to rank peaks
res <- list()

for(j in 7:length(patients)){
patient <- patients[j]
# get nf cutsite files
nf_cutsites <- list.files(path = "cutsites",
                          pattern = paste0("EPICC_",patient,".*_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz$"),
                          full.names = TRUE,
                          recursive=TRUE)
nf_cutsites <- nf_cutsites[!grepl("_E",nf_cutsites)]
  
# take mutations which are subclonal for this patient 
my_muts <- muts[mcols(muts)[,patient]=="subclonal"]

for(i in 1:length(my_muts)){
  mut <- my_muts[i]
  # get enhancer region
  ehs <- mut
  start(ehs) <-  as.integer(gsub("-.*","",ehs$Intervals))
  end(ehs) <-  as.integer(gsub(".*-","",ehs$Intervals))
  ehs <- makeGRangesFromDataFrame(ehs,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  starts.in.df.are.0based=FALSE)
  
  # get mutated samples
  mutated <- paste0(gsub(";","_|",mut$MutAllSam),"_")
  mutated_files <- nf_cutsites[grepl(mutated,nf_cutsites)]
  
  # get WTs taking only those with WGS
  wt_samples <- matches[matches$Patient==patient & !(matches$Sample %in% unlist(strsplit(mut$MutAllSam,split=";"))) ,"Sample"]
  wt_files <- nf_cutsites[grepl(paste0(paste0(wt_samples,collapse="_|"),"_") ,nf_cutsites)]
  
  # import the cut sites as GRanges
  mutated_gr <- unlist(lapply(1:length(mutated_files), function(s) {
    gr <-import(mutated_files[s],
                format="BED",
                genome="hg38",
                which=ehs) 
    len <- length(gr)
    #normalise the length by reads per million in peaks 
    sample <- gsub("_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz","",basename(mutated_files[s]))
    mil_mut <- sum(counts[row.names(peak_matches_df)[peak_matches_df[,patient]],sample])/1000000
    len <-  len/mil_mut
  }))

  wt_gr <-  unlist(lapply(1:length(wt_files), function(s) {
    gr <-import(wt_files[s],
                format="BED",
                genome="hg38",
                which=ehs) 
    len <- length(gr)
    #normalise the length by reads per million in peaks
    sample <- gsub("_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz","",basename(wt_files[s]))
    mil_wt <- sum(counts[row.names(peak_matches_df)[peak_matches_df[,patient]],sample])/1000000
    len <-  len/mil_wt
  }))
  
  # perform T-test
  test <- NA
  try(test <- t.test(mutated_gr, wt_gr, alternative = "two.sided", var.equal = FALSE))
  try(test <- test$p.value)
  
  # get total read numbers for the enhancer
  nmut <- length(mutated_files)
  nwt <- length(wt_files)
  # get mean read numbers and logFC for the enhancer
  readsmut <- mean(mutated_gr)
  readswt <- mean(wt_gr)
  logFC <- log2(readsmut/readswt)

  # get total read numbers for the mutated/wt samples
  total_reads_mut <- sum(counts[,gsub("_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz","",basename(mutated_files))])
  total_reads_wt <- sum(counts[,gsub("_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz","",basename(wt_files))])
  # get total read numbers within peaks for the mutated/wt samples
  total_readsip_mut <- sum(counts[row.names(peak_matches_df)[peak_matches_df[,patient]],
                                  gsub("_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz","",basename(mutated_files))])
  total_readsip_wt <- sum(counts[row.names(peak_matches_df)[peak_matches_df[,patient]],
                                 gsub("_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz","",basename(wt_files))])
  
  # put the metrics together to form results table
  res[[length(res)+1]] <- c(mut$Locus,mut$Gene,mut$Ensembl,patient,test,logFC,mut$Mutes,
                            nmut,nwt,readsmut,readswt,total_reads_mut,total_reads_wt,
                            total_readsip_mut,total_readsip_wt,
                            mut$HWrec,mut$HWRNArec,mut$MutAllSam,wt_samples
                            )
   
}
 
}


res.df <- do.call("rbind",res)
colnames(res.df) <- c("ID","Gene","Ensembl","patient","p.value_one-sided","LogFC","EffectSize",
                      "n_mut","n_wt","CPM_mut_mean","CPM_wt_mean","total_reads_mut","total_reads_wt",
                      "total_reads_in_peaks_mut","total_reads_in_peaks_wt",
                      "HWrec","HWRNArec","mut_samples","wt_samples")
write.table(res.df,"t_test.csv",sep="\t", quote=F, row.names=F)

