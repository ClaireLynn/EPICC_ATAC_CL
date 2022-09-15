library(rtracklayer)
library(edgeR)
library(GenomicFeatures)
library(GenomicRanges)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(biomaRt)
library(Gviz)
library(caTools)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(scales)
library(gtable)

library(lattice)

source("horizon_plot_functions_commonE_CL6.R")
results <- readRDS("all_deseq_results.rds") 

cov <- readRDS("../recurrent_peaks_final/cpm_peaknorm_cov_nf_1bp.rds")
beds <- list.files(path = "~/Dropbox (ICR)/EPICC/EPICC/prospective_cohort/data_analysis/ATAC_peak_calls/",
                   pattern = ".*bed$",
                   full.names = TRUE,
                   recursive=TRUE)
beds <- beds[grepl("_A|_B|_C|_D|C542_nucleosome_free-region_F", beds)]

beds_gr <- lapply(1:length(beds), function(bed){
  bed_gr <- makeGRangesFromDataFrame(read.table(beds[bed], header=T),
                                    keep.extra.columns=TRUE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     starts.in.df.are.0based=FALSE)
  bed_gr <- subsetByOverlaps(bed_gr,peaks)
  return(bed_gr)
})
names(beds_gr) <- apply(str_match(basename(beds),"(.*?)_.*region_(.*?)-GRCh38.*")[,c(2,3)],
                        1 , paste , collapse = "_" )

# For each recurrent peak / patient combo, are there any peaks,
# if not DON'T plot
peak_matches <- lapply(1:length(beds), function(bed){
  peak_match <- makeGRangesFromDataFrame(read.table(beds[bed], header=T),
                                     keep.extra.columns=TRUE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     starts.in.df.are.0based=FALSE)
  peak_match <- subsetByOverlaps(peaks, peak_match)
  peak_match <-mcols(peak_match)$peak
  return(peak_match)
})
names(peak_matches) <- apply(str_match(basename(beds),"(.*?)_.*region_(.*?)-GRCh38.*")[,c(2,3)],
                        1 , paste , collapse = "_" )


peak_matches_df <- matrix(nrow=length(peaks$peak) , ncol= length(patients), FALSE)
rownames(peak_matches_df) <- peaks$peak
colnames(peak_matches_df) <- patients
for(patient in patients){
  peaks_with_peaks<-mcols(peaks[peaks$peak %in% unique(unlist(peak_matches[grepl(patient,names(peak_matches))])),])$peak
  peak_matches_df[peaks_with_peaks,patient] <- TRUE
}



#Get peaks again for tracks
beds <- list.files(path = "~/Dropbox (ICR)/EPICC/ATAC_CL/peaks",
                   pattern = ".*bed$",
                   full.names = TRUE,
                   recursive=TRUE)
beds <- beds[grepl("_A|_B|_C|_D|_E|C542_nucleosome_free-region_F", beds)]

beds_gr <- lapply(1:length(beds), function(bed){
  bed_gr <- makeGRangesFromDataFrame(read.table(beds[bed], header=T),
                                     keep.extra.columns=TRUE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     starts.in.df.are.0based=FALSE)
  bed_gr <- subsetByOverlaps(bed_gr,peaks)
  return(bed_gr)
})
names(beds_gr) <- apply(str_match(basename(beds),"(.*?)_.*region_(.*?)-GRCh38.*")[,c(2,3)],
                        1 , paste , collapse = "_" )

peak_tracks <- lapply(1:length(beds_gr),function(i){
  region=unlist(strsplit(names(beds_gr[i]),split="_"))[2]
  color <- col[region]
  return(AnnotationTrack(beds_gr[[i]],name=NULL,fill=color))})
names(peak_tracks) <- names(beds_gr)


saveRDS(peak_matches_df,"~/subclonal_redo/peak_matches_df.rds")

#MAKE SUBCLONAL MATRIX FOR HEATMAP
sc_res <- results[,grep("region_nf_padj",colnames(results))]
sc_res <- results[,paste0(patients,"_region_nf_padj")]
colnames(sc_res) <- patients
sc_res<-mapply(sc_res, FUN=as.numeric)
rownames(sc_res) <- results$peak
#get the fcs
fcc <- fc
fcc[is.na(fcc)] <- 0
fcc<-fcc[,grepl(".pure",colnames(fcc))]
colnames(fcc) <- unlist(lapply(strsplit(colnames(fcc),"\\."),"[",1))
fcc<-fcc[,colnames(sc_res)]
fccl <- fcc < 0
fccg <- fcc > 0
gains <- results[results$event_type_chr=="gain","peak"]
losses <- results[results$event_type_chr=="loss","peak"]
subclonal_matrix <- sc_res < 0.05
subclonal_matrix[is.na(subclonal_matrix)] <- FALSE
# If it is a loss the direction must be the same!
subclonal_matrix[losses,] <-(subclonal_matrix[losses,] + fccl[losses,]) >1
# If it is a gain, there must be at least 1 peak & it must be in the right direction!
subclonal_matrix[gains,] <- (subclonal_matrix[gains,] + peak_matches_df[gains,] + fccg[gains,]) > 2

cancer_s_mat <-cancer_s_matrix[peaks$peak,]
colnames(cancer_s_mat) <-  gsub("\\.pure","",colnames(cancer_s_mat))
cancer_s_mat <-cancer_s_mat[,colnames(subclonal_matrix)]

both <- cancer_s_mat+subclonal_matrix


saveRDS(subclonal_matrix,"~/git_projects/EPICC_ATAC/subclonal_matrix.rds")

library(ggplot2)

# Plot subclonal
for(i in 1:length(patients)){
  peakss <- peaks
  patient <- patients[i]
  peaks_with_peaks<-mcols(peaks[peaks$peak %in% unique(unlist(peak_matches[grepl(patient,names(peak_matches))])),])$peak
  for(t in 1:length(peakss)){
    peak <- peakss$peak[t]
    present <-ifelse(peak %in% peaks_with_peaks, "gotpeaks","nopeaks")
    logfc <-round(fc[peak,paste0(patient,".pure")],digits = 2)
    sig <- sig_filter[peak,paste0(patient,".pure")]
    status<-get_status2(patient=patient,
                        peak=peakss$peak[t])
    try(plot_info(pdf_name=paste0("~/subclonal_redo/plots2/",patients[i],"_",
                                  peakss$id[t],"_",
                                  peakss$type[t],"_",
                                  peakss$event_type_chr[t],"_",
                                  status,"_",logfc,"_",sig,"_",present,"_line.pdf"),
                  patient=patient,
                  peak=peakss$peak[t],
                  type="l"))
    graphics.off()
  }}
