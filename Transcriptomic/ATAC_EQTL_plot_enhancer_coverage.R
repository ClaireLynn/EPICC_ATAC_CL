# This script gets the coverage of enhancers for specified samples (WT/mutant),
# performs normalisation and plots the tracks

library(rtracklayer)
library(RColorBrewer)
library(Gviz)

# set smoothing and order of tracks
args = commandArgs(trailingOnly=TRUE)
SMOOTH <- as.numeric(args[1])
MAX <- as.character(args[2])
dir.create(paste0("plots_",MAX,"_SMOOTH_",SMOOTH))

# set colours
col <- brewer.pal(name="Set1", n=9)[c(1,2,5)]
names(col) <- c("mutated","wt","normal")

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

# compare samples from the same case with/without the mutation
for(j in 1:length(patients)){
patient <- patients[j]
# get nf cutsite files
nf_cutsites <- list.files(path = "cutsites",
                          pattern = paste0("EPICC_",patient,".*_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz$"),
                          full.names = TRUE,
                          recursive=TRUE)
nf_cutsites <- nf_cutsites[!grepl("_E",nf_cutsites)]

# take mutations which are subclonal for this patient 
my_muts <- muts[mcols(muts)[,patient]=="subclonal"]

# for each mutation, plot WT and mutant ATAC
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
  ehs_resize <-resize(ehs, width = width(ehs)+(5000*2), fix = "center")
  
  # get 1bp bins across enhancers so that coverage has 1bp resolution
  bins <- unlist(tile(ehs_resize,width=1))
  
  # get mutated samples
  mutated <- paste0(gsub(";","_|",mut$MutAllSam),"_")
  mutated_files <- nf_cutsites[grepl(mutated,nf_cutsites)]

  # get WTs taking only those with WGS
  wt_samples <- matches[matches$Patient==patient & !(matches$Sample %in% unlist(strsplit(mut$MutAllSam,split=";"))) ,"Sample"]
  wt_files <- nf_cutsites[grepl(paste0(paste0(wt_samples,collapse="_|"),"_") ,nf_cutsites)]
  
  # import the cut sites as GRanges and extend them to nucleosome size
  mutated_gr <- unlist(GRangesList(lapply(1:length(mutated_files), function(s) {
    gr <-import(mutated_files[s],
                format="BED",
                genome="hg38",
                which=ehs_resize) 
    gr <- resize(gr, width = 100, fix = "center")
    return(gr)
  })))
  
  wt_gr <-  unlist(GRangesList(lapply(1:length(wt_files), function(s) {
    gr <-import(wt_files[s],
                format="BED",
                genome="hg38",
                which=ehs_resize) 
    gr <- resize(gr, width = 100, fix = "center")
    return(gr)
  })))
  
  # get cutsites for pooled normal
  gr <-import("cutsites/region_E_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz",
                format="BED",
                genome="hg38",
                which=ehs_resize) 
  normal_gr <- resize(gr, width = 100, fix = "center")
  # normalise coverage per million reads within peaks
  # 43732335 is total reads in peaks from megabulk_data.rds
  millions_E <- 43732335/1000000
  
  # get reads in peaks for mutant/wt samples
  mutants <- unlist(strsplit(mut$MutAllSam, split=";"))
  mutants <- mutants[grepl(patient,mutants)]
  mutcounts <- counts[,grepl(paste0(mutants,collapse="|"),colnames(counts)),drop=FALSE]
  wtcounts <- counts[,grepl(paste0(wt_samples,collapse="|"),colnames(counts)),drop=FALSE]
  # get peak matches
  match <- setNames(peak_matches_df[,patient],rownames(peak_matches_df))
  match <- names(match[which(match==TRUE)])
  mutcounts <- mutcounts[match,]
  wtcounts <- wtcounts[match,]
  # normalise coverage per million reads within peaks
  libsize <- c(sum(mutcounts),sum(wtcounts))
  names(libsize) <- c("mut","wt")
  millions <- libsize/1000000
  mutated_cov <- coverage(mutated_gr)/millions["mut"]
  mutated_cov <- mutated_cov[seqlevels(bins)]
  mutated_cov <- GenomicRanges::binnedAverage(bins, mutated_cov, "score")
  wt_cov <-  coverage(wt_gr)/millions["wt"]
  wt_cov <- wt_cov[seqlevels(bins)]
  wt_cov <- GenomicRanges::binnedAverage(bins, wt_cov, "score")
  normal_cov <-  coverage(normal_gr)/millions_E
  normal_cov <- normal_cov[seqlevels(bins)]
  normal_cov <- GenomicRanges::binnedAverage(bins, normal_cov, "score")

  # get the order of presentation so that max value is plotted first, aesthetic only
  if(MAX=="peak"){
  max<-c(max(mutated_cov$score),max(wt_cov$score),max(normal_cov$score))
  }
  if(MAX=="sum"){
  max<-c(sum(mutated_cov$score),sum(wt_cov$score),sum(normal_cov$score))
  }

  max[which(!is.finite(max))] <- 0
  names(max) <- c("mutated_cov","wt_cov","normal_cov")
  max_order <- order(max, decreasing=T)
  max_order <- names(max[max_order])

  # get tracks from coverage and colour them
  mutated_track <- DataTrack(range=mutated_cov,
                             col.line=col,
                             from=start(ehs_resize),
                             to = end(ehs_resize),
                             chr=as.character(seqnames(ehs_resize))
                                       )
  displayPars(mutated_track) <- list(groups=factor("mutated", 
                                              levels = c("mutated","wt","normal")),
                                col=col["mutated"],
                                legend = TRUE,
                                missingAsZero=F,
                                lwd=3)
  
  wt_track <- DataTrack(range=wt_cov,
                        col.line=col,
                        from=start(ehs_resize),
                        to = end(ehs_resize),
                        chr=as.character(seqnames(ehs_resize)) 
                             )
  displayPars(wt_track) <- list(groups=factor("wt", 
                                           levels = c("mutated","wt","normal")),
                             col=col["wt"],
                             legend = TRUE,
                             missingAsZero=F,
                             lwd=3)
  
  normal_track <- DataTrack(range=normal_cov,
                        col.line=col,
                        from=start(ehs_resize),
                        to = end(ehs_resize),
                        chr=as.character(seqnames(ehs_resize)) 
  )
  displayPars(normal_track) <- list(groups=factor("normal", 
                                              levels = c("mutated","wt","normal")),
                                col=col["normal"],
                                legend = TRUE,
                                missingAsZero=F,
                                lwd=3)
 
  # overlay tracks for plotting
  tracks <-list(mutated_track,wt_track,normal_track)
  names(tracks)=c("mutated_cov","wt_cov","normal_cov")
  tracks <- tracks[max_order]
  names(tracks[[1]]) <- "CPM"
  ylims <- range(lapply(tracks,values))
  tracks_to_plot<-OverlayTrack(tracks,name=as.character(patient))

 locus <- gsub("/","_",mut$Locus)
 
 # Add highlight for the location of the mutation within the enhancer
 tracks_to_plot_hl <- HighlightTrack(trackList = list(tracks_to_plot),
                     range=my_muts[i],col="#E41A1C")

  pdf(paste0("plots_",MAX,"_SMOOTH_",SMOOTH,"/",patient,"_",locus,"_",mut$Gene,".pdf"),height=4.5,width=11)
  plotTracks(background.title="white",
             fontcolor.title="black",
             col.axis="black",
             showTitle=T, clip=FALSE,
             c(GenomeAxisTrack(labelPos = "above",name=NULL,add53 = F, add35 = F),
               AnnotationTrack(ehs, name=NULL,fill="grey", group=paste(mut$Gene,"enhancer"),just.group="above",groupAnnotation="group"),
               AnnotationTrack(my_muts[i], name = NULL, fill="#E41A1C",col="#E41A1C",
                               group=paste0("  ",gsub("/",">",gsub("_"," ",mut$Locus))),groupAnnotation="group",
                               fontcolor.group = "#e65052",just.group = "right"),
                               tracks_to_plot_hl),
             legend=T,
             type="l",
             from=start(ehs_resize),
             to = end(ehs_resize),
             min.distance=10,
             chromosome=as.character(seqnames(ehs_resize)),
             lwd.baseline=1,
             col.baseline="black",
             windowSize =SMOOTH,window=-1,
             ylim=ylims,cex.title=1,cex.axis=0.7, cex.main=1, add=TRUE,
             sizes=c(0.2,0.1,0.1,1),
             detail="coverage")
dev.off()

  }

}








