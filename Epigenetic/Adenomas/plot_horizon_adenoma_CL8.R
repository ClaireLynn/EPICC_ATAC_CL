peaks$symbol <- str_replace_all(peaks$symbol, "−", "-")

beds <- list.files(path = "~/Dropbox (ICR)/EPICC/ATAC_CL/peaks",
                   pattern = ".*bed$",
                   full.names = TRUE,
                   recursive=TRUE)

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

saveRDS(peak_matches,"peak_matches_adenoma.rds")

peak_tracks <- lapply(1:length(beds_gr),function(i){
  region=unlist(strsplit(names(beds_gr[i]),split="_"))[2]
  color <- col[region]
  return(AnnotationTrack(beds_gr[[i]],name=NULL,fill=color))})
names(peak_tracks) <- names(beds_gr)

saveRDS(peak_matches,"peak_matches_adenoma.RDS")
saveRDS(beds_gr,"beds_gr_adenoma.RDS")
saveRDS(peak_tracks,"peak_tracks_adenoma.RDS")


sig_filter <- final_reccurence_count$sig_matrix
peaks$symbol <- str_replace_all(peaks$symbol, "−", "-")
for(i in 1){
  patient <- as.character(patients[i])
  colname<-paste0(patient,".adenoma")
  peakss <- peaks
  peaks_with_peaks<-mcols(peaks[peaks$peak %in% unique(unlist(peak_matches[grepl(patient,names(peak_matches))])),])$peak
  for(t in 1:length(peaks)){
    peak <- peakss$peak[t]
    present <-ifelse(peak %in% peaks_with_peaks, "gotpeaks","nopeaks")
    NvA_logfc <- round(fc[peak,colname], digits=2)
    sigNvA <- sig_filter[peak,colname]
    status<-get_status_adenoma(patient=patient,
                       peak=peakss$peak[t])
    try(plot_info(pdf_name=paste0("~/subclonal_redo/adenoma/",patient,"_",
                                  peakss$id[t],"_",
                                  peakss$type[t],"_",
                                  peakss$event_type_chr[t],"_",
                                  status,"_NvA.",NvA_logfc,"_",sigNvA,"_",
                                  present,"_line.pdf"),
                  patient=patient,
                  peak=peakss$peak[t],
                  type="l"))
    graphics.off()
  }
  
  
  
  }


for(i in 7){
  patient <- as.character(patients[i])
  peakss <- peaks
  peaks_with_peaks<-mcols(peaks[peaks$peak %in% unique(unlist(peak_matches[grepl(patient,names(peak_matches))])),])$peak
  for(t in 1:length(peaks)){
    peak <- peakss$peak[t]
    present <-ifelse(peak %in% peaks_with_peaks, "gotpeaks","nopeaks")
    status<-get_status_adenoma_C561(peak=peakss$peak[t])
    try(plot_info_C561(pdf_name=paste0("~/subclonal_redo/adenoma2/",patient,"_",
                                  peakss$id[t],"_",
                                  peakss$type[t],"_",
                                  peakss$event_type_chr[t],"_",
                                  status,"_",
                                  present,"_line.pdf"),
                  patient=patient,
                  peak=peakss$peak[t],
                  type="l"))
    graphics.off()
  }
  
  
  
}




