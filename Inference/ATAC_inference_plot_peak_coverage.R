library(rtracklayer)
library(Gviz)
library(RColorBrewer)
library(tidyverse)

# Pval table function
#######################################################

Pval_table <- function(res,peak=peak){
  results_plot <- res[res$peak==peak,]
  results_table<- data.frame(logFC=signif(as.numeric(results_plot[,grepl("log2FoldChange",colnames(results_plot))]),digits=2),
                             padj=signif(as.numeric(results_plot[,grepl("padj",colnames(results_plot))]),digits=2))
  row.names(results_table) <- unique(gsub(".*\\.","",grep("log2FoldChange",colnames(results_plot),value=T)))
  return(results_table)
}

# Track plotting function
#######################################################

plot_horizon <- function(cov,col,patient,peak,smooth,lwd=4,type="l",legend){
  peak_gr <- resizePeaks[resizePeaks$peak == peak,]
  start=start(resizePeaks[resizePeaks$peak == peak,])
  end=end(resizePeaks[resizePeaks$peak == peak,])
  chr=as.character(as.data.frame(resizePeaks)[resizePeaks$peak == peak,1])
  
#Get the order of presentation so that max value is plotted first
  cov<-lapply(cov,subsetByOverlaps,peak_gr)
  suppressWarnings(max<-unlist(lapply(1:length(groups),function(n){
    c <- groups[n]
    max(cov[[c]]$score)})))
  
  max[which(!is.finite(max))] <- 0
  names(max) <- groups
  max_order <- order(max, decreasing=T)
  max_order <- names(max[max_order])

# get tracks from coverage and colour them
  tracks <- lapply(1:length(groups), function(c) {
    mycov <- cov[[c]]
    group <- groups[c]
    track <- DataTrack(range=mycov,
                       fill.mountain="transparent",
                       col.mountain=col,
                       name=groups[c],fill=col,
                       col.histogram=col,
                       fill.histogram="transparent",
                       fill="transparent",
                       col=col,
                       col.line=col,
                       from=start,
                       to = end,
                       chr=chr
                       )
    displayPars(track) <- list(groups = factor(group,
                                               levels = groups),
                               col=col[group],
                               legend = TRUE,
                               missingAsZero=F,
                               lwd=3)
    return(track)})
  
  names(tracks) <- groups
  tracks <- tracks[max_order]
  names(tracks[[1]]) <- paste0(as.character(patient)," normalised coverage")

  ylims <- range(lapply(tracks,values))
  tracks_to_plot<-OverlayTrack(tracks,name=as.character(patient),
                               ylim=ylims)
  
  plotTracks(background.title="white",fontcolor.title="black",col.axis="black",
             showTitle=T, clip=FALSE,
             c(AnnotationTrack(peaks, name = NULL, fill="grey"),
               tracks_to_plot),
             type="l",
             from=start,
             to = end,
             min.distance=10,
             chromosome=chr,
             lwd.baseline=1,
             col.baseline="black",
             windowSize =40,window=-1,
             main=paste(mcols(peaks)[mcols(peaks)$peak==peak,"id"],
                        mcols(peaks)[mcols(peaks)$peak==peak,"type"],"Padj:",
                        round(mcols(peaks)[mcols(peaks)$peak==peak,"padj"],digits=3),"LogFC:",
                        round(mcols(peaks)[mcols(peaks)$peak==peak,"log2FoldChange"],digits=3),"\n",
                        peak,
                        sep=" "),
             ylim=ylims,cex.title=1,cex.axis=0.7, cex.main=1, add=TRUE,
             sizes=c(0.02,1),
             detail="coverage")
}

# Plot significant peaks for all 
#######################################################

patients <- c("C525","C538","C539","C542","C549","C516","C518","C524","C531","C551","C559","C562")

for(pat in 1:length(patients)){
  patient <- patients[pat]
  coldata <- readRDS("coldata_inf.rds")
# Get patient and remove samples with no purity info/ adenomas from coldata
  coldata <- coldata[coldata$patient==patient,]
  coldata <- coldata[!(is.na(coldata$purity) | is.na(coldata$group) | coldata$purity == 0 | coldata$subtissue == "adenoma"), ]
  multi <- length(unique(coldata$group)) > 2
  res_files <- list.files(path = "fresh_effect_of_group_nfreads_LRT", pattern=paste0(patient,".*.rds"),full.names =TRUE)
  
# If there are multiple comparisons
  if(multi==TRUE){
    res <- lapply(res_files,readRDS)
    for(x in 1:length(res_files)){
      colnames(res[[x]])[9:15] <- paste0(colnames(res[[x]])[9:15],
                                         ".", gsub("_.*","",basename(res_files[x])))
    }
    for(x in 2:length(res_files)){
      res[[x]] <- res[[x]][,c(1,9:15)]
    }
    res <- res %>% reduce(left_join, by = "peak")
  }else{
    res <- readRDS(res_files)
  }
  
  chr <- unlist(lapply(str_split(res[res$padj<0.1,"peak"], ":"), `[[`, 1))
  start <- unlist(lapply(str_split(lapply(str_split(res[res$padj<0.1,"peak"], ":"), `[[`, 2), "-"),  `[[`, 1))
  end <- unlist(lapply(str_split(lapply(str_split(res[res$padj<0.1,"peak"], ":"), `[[`, 2), "-"),  `[[`, 2))
  df <- as.data.frame(cbind(chr,start,end,res[res$padj<0.1,]))
  row.names(df) <- NULL
  
  summary<-readRDS("summary.rds")
  df <- merge(summary[,c(1,2)],df[,-5],by.x=1,by.y=4)
  peaks <- makeGRangesFromDataFrame(df,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    starts.in.df.are.0based=FALSE)
  resizePeaks <- resize(peaks, width = 2000, fix = "center")
  nf_cutsites <- list.files(path = "../cutsites",
                            pattern = paste0("EPICC_",patient,".*_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz$"),
                            full.names = TRUE,
                            recursive=TRUE)
  
  samples <- str_match(nf_cutsites, "EPICC_\\s*(.*?)\\s*_C1_nuc")[,2]
  include <- samples %in% rownames(coldata)
  nf_cutsites <- nf_cutsites[include]
  groups<-levels(as.factor(coldata$group))
  
  sample <- str_match(nf_cutsites, "EPICC_\\s*(.*?)\\s*_C1_nuc")[,2]
  files <- as.data.frame(cbind(nf_cutsites,sample))
  files$group <- coldata[files$sample,"group"]
  
  invisible(lapply(1:length(groups), function(c){
    gr <- unlist(GRangesList(
      lapply(which(files$group==groups[c]), function(s) {
        gr <-import(nf_cutsites[s],
                    format="BED",
                    genome="hg38",
                    which=resizePeaks)
        gr <- resize(gr, width = 100, fix = "center")
        return(gr)
        })))
    saveRDS(gr,paste0("peak_plots/",groups[c],"_sigpeaks_reads.rds"))
    }))
  
  bins <- unlist(tile(resizePeaks,width=1))

# Get normalisation
  counts2 <- "atacseq_counts_per_peak.rds"
  counts2 <- readRDS(counts2)
  colnames(counts2) <- gsub("EPICC_|_C1$","",colnames(counts2))
  counts2 <- counts2[,rownames(coldata)]
# get peak matches
  peak_matches_df<-readRDS("peak_matches_df.rds")
  match <- setNames(peak_matches_df[,patient],rownames(peak_matches_df))
  match <- names(match[which(match==TRUE)])
  counts2 <- counts2[match,]
  
  libsize <- unlist(lapply(1:length(groups), function(c){
    samples<-row.names(coldata[coldata$group==groups[c],])
    sum(counts2[,samples])
  }))
  names(libsize) <- groups

  # write the normalised coverage to a file
  cov <- function(libsize,name){
    millions <- libsize / 1000000
    covs <- lapply(1:length(groups), function(c) {
      cov <- coverage(readRDS(paste0("peak_plots/",groups[c],"_sigpeaks_reads.rds")))/millions[c]
      cov <- cov[seqlevels(bins)]
      cov <- GenomicRanges::binnedAverage(bins, cov, "score")
      return(cov)})
    names(covs) <- groups
    saveRDS(covs,name)
    }

  cov(libsize,paste0("peak_plots/",patient,"_cpm_sigpeaks_reads_cov.rds"))

# get normal
  normal_gr <-import("../cutsites/region_E_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz",
                     format="BED",
                     genome="hg38",
                     which=resizePeaks)
  normal_gr <- resize(normal_gr, width = 100, fix = "center")
  millions_E <- 43732335/1000000
  normal_cov <-  coverage(normal_gr)/millions_E
  normal_cov <- normal_cov[seqlevels(bins)]
  normal_cov <- GenomicRanges::binnedAverage(bins, normal_cov, "score")
  
  col <- brewer.pal(name="Set1", n=9)[c(1:length(groups),5)]
  names(col) <- c(groups,"normal")

  groups <- c(groups,"normal")

  cov <- readRDS(paste0("peak_plots/",patient,"_cpm_sigpeaks_reads_cov.rds"))
  cov <- c(cov,normal_cov)
  names(cov) <- groups


  dir.create(paste0("peak_plots/",patient))
  peak_name <-  gsub("\\)","",gsub(", | \\(","_", peaks$id))
  
# Plot subclonal
# If there is more than one comparison, make a table with pvalue and logFC
  if(multi==TRUE){
    for(t in 1:length(peaks$peak)){
      peak <- peaks$peak[t]
      pdf(paste0("peak_plots/",patient,"/",patient,"_",
                 peak_name[t],"_lineE.pdf"), height=4, width=6)
      vp <- viewport(x = 0.87, y = 0.9,
                     width = 0.2, height = 0.2,
                     just = c("right", "top"))
      pushViewport(vp)
      grid.draw(tableGrob(Pval_table(res,peak=peak)))
      popViewport()
      print(plot_horizon(cov=cov,
                         col=col,
                         patient=patient,
                         peak=peak,
                         smooth=30,
                         lwd=1.5,
                         legend=T,
                         type="l"))
      dev.off()
    }} else{
  for(t in 1:length(peaks$peak)){
    peak <- peaks$peak[t]
    pdf(paste0("peak_plots/",patient,"/",patient,"_",
               peak_name[t],"_",round(peaks$padj[t],digits=3),"_lineE.pdf"), height=4, width=6)
    print(plot_horizon(cov=cov,
                       col=col,
                       patient=patient,
                       peak=peak,
                       smooth=30,
                       lwd=1.5,
                       legend=T,
                       type="l"))
    dev.off()
  }}
  }


